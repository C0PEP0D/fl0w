#include "getjhtdb/getData.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <curl/curl.h>
#include "cjson/cJSON.h"
#include <math.h>

char *convert_points_to_string(float points[][3], int num_points) {
	// Estimate initial buffer size based on typical point size
	size_t initial_size = num_points * 64; // 64 is an inital estimate...
	char *points_str = (char *)malloc(initial_size);
	if (points_str == NULL) {
		fprintf(stderr, "Failed to allocate memory for points_str\n");
		return NULL;
	}

	// Initialize empty string
	points_str[0] = '\0';

	// Keep track of the current buffer size and used size
	size_t current_size = initial_size;
	size_t used_size = 0;

	for (int i = 0; i < num_points; i++) {
		char coord_str[64];
		int len = snprintf(coord_str, sizeof(coord_str), "%.8f\t%.8f\t%.8f\n", points[i][0], points[i][1], points[i][2]);

		// Check if there's enough space left in the buffer
		if (used_size + len + 1 > current_size) { // +1 for null terminator...
			// Double the buffer size
			current_size *= 2;
			char *new_points_str = realloc(points_str, current_size);
			if (new_points_str == NULL) {
				fprintf(stderr, "Failed to reallocate memory for points_str\n");
				free(points_str);
				return NULL;
			}
			points_str = new_points_str;
		}

		// Append the coordinate string to the buffer
		strcat(points_str, coord_str);
		used_size += len;
	}

	return points_str;
}

size_t write_callback(void *contents, size_t size, size_t nmemb, void *userp) {
	size_t real_size = size * nmemb;
	char **response_ptr = (char **)userp;

	// Reallocate the memory with the new size
	*response_ptr = realloc(*response_ptr, strlen(*response_ptr) + real_size + 1);

	if (*response_ptr == NULL) {
		// Out of memory, error message
		fprintf(stderr, "Not enough memory (realloc returned NULL)\n");
		return 0;
	}

	// Append the new data to the response buffer
	strncat(*response_ptr, (char *)contents, real_size);

	return real_size;
}

double **getData(const char *auth_token, const char *dataset, float time, const char *spatial_interpolation, 
				 const char *temporal_interpolation, const char *variable, const char *spatial_operator, 
				 float points[][3], int num_points) {
	CURL *curl;
	CURLcode res;
	
	char *points_str = convert_points_to_string(points, num_points);

	// Allocate initial memory for the response
	char *response = (char *)malloc(1);

	// EMPTY STRING
	response[0] = '\0';  
	
	char url[4096];
	int num_columns;

	char time_str[50];
	snprintf(time_str, sizeof(time_str), "%.6f", time);
	
	if (strcmp(variable, "velocity") == 0 || strcmp(variable, "magneticfield") == 0 || strcmp(variable, "vectorpotential") == 0 || strcmp(variable, "force") == 0 || strcmp(variable, "position") == 0) {
		if (strcmp(spatial_operator, "field") == 0) {
	  num_columns = 3;
		
		} else if (strcmp(spatial_operator, "gradient") == 0) {
			num_columns = 9;
		} else if (strcmp(spatial_operator, "hessian") == 0) {
			num_columns = 18;
		}
	} else if (strcmp(variable, "pressure") == 0 || strcmp(variable, "density") == 0 || strcmp(variable, "energy") == 0 || strcmp(variable, "temperature") == 0) {
		if (strcmp(spatial_operator, "field") == 0) {
			num_columns = 1;
		} else if (strcmp(spatial_operator, "gradient") == 0) {
			num_columns = 3;
		} else if (strcmp(spatial_operator, "hessian") == 0) {
			num_columns = 6;
		}
	}

	curl_global_init(CURL_GLOBAL_DEFAULT);
	curl = curl_easy_init();

	double **data_operator_2D = NULL;    

	//printf("check point none and none %d %s %s", num_points, time_end, delta_t);
	  
	data_operator_2D = (double **)malloc(num_points * sizeof(double *));
	
	for (int i = 0; i < num_points; i++) {
	  data_operator_2D[i] = (double *)malloc(num_columns * sizeof(double));
	}
	
	if (curl) {
		// Construct the URL
		snprintf(
			url, 
			sizeof(url), 
			"https://web.idies.jhu.edu/turbulence-svc-test/values?authToken=%s&dataset=%s&function=GetVariable&var=%s&t=%s&sint=%s&sop=%s&tint=%s",
			auth_token, 
			dataset, 
			variable, 
			time_str, 
			spatial_interpolation, 
			spatial_operator, 
			temporal_interpolation
		);
		
		// Set HTTP request options
		curl_easy_setopt(curl, CURLOPT_URL, url);
		curl_easy_setopt(curl, CURLOPT_POSTFIELDS, points_str);
		curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_callback);
		curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);

		res = curl_easy_perform(curl);

	  
		if (res != CURLE_OK) {
			fprintf(stderr, "curl_easy_perform() failed: %s\n", curl_easy_strerror(res));
		} else {
			// Process the JSON response
			cJSON *json = cJSON_Parse(response);
	
			// Check if there is an error in the response
			cJSON *title = cJSON_GetObjectItem(json, "title");
			if (title && strcmp(title->valuestring, "Error") == 0) {
				cJSON *description = cJSON_GetObjectItem(json, "description");
	  			if (description) {
					printf("Error is : %s\n", description->valuestring);
	  			}
			} else {
	  			int size = cJSON_GetArraySize(json);
	  			
				  for (int p = 0; p < size; p++) {
					cJSON *json_elem = cJSON_GetArrayItem(json, p);
					for (int j = 0; j < num_columns; j++) {
						data_operator_2D[p][j] = cJSON_GetNumberValue(cJSON_GetArrayItem(json_elem, j));
					}
				}
			}

			cJSON_Delete(json);
		}
		
		curl_easy_cleanup(curl);
	}
	
	free(points_str);
	curl_global_cleanup();

	return data_operator_2D;
}

void writeVTK(const char *filename, double **data, int nx, int ny, int nz, int num_columns, float points[][3]) {
	FILE *fp = fopen(filename, "w");
	if (!fp) {
		fprintf(stderr, "Cannot open file %s for writing\n", filename);
		return;
	}

	// VTK Header, lots of stuff need to be added here.....
	fprintf(fp, "# vtk DataFile Version 3.0\n");
	fprintf(fp, "VTK output\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET STRUCTURED_GRID\n");
	fprintf(fp, "DIMENSIONS %d %d %d\n", nx, ny, nz);
	fprintf(fp, "POINTS %d float\n", nx * ny * nz);

	// Write the points coordinates
	for (int i = 0; i < nx * ny * nz; i++) {
		fprintf(fp, "%f %f %f\n", points[i][0], points[i][1], points[i][2]);
	}

	// Write the data
	fprintf(fp, "POINT_DATA %d\n", nx * ny * nz);
	for (int col = 0; col < num_columns; col++) {
		fprintf(fp, "SCALARS data_%d float 1\n", col);
		fprintf(fp, "LOOKUP_TABLE default\n");
		for (int i = 0; i < nx * ny * nz; i++) {
			fprintf(fp, "%f\n", data[i][col]);
		}
	}

	fclose(fp);
}
