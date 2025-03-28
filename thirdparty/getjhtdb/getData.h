#ifndef GETDATA_H
#define GETDATA_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <curl/curl.h>
#include "cjson/cJSON.h"

char *convert_points_to_string(float points[][3], int num_points);
size_t write_callback(void *contents, size_t size, size_t nmemb, void *userp);
double **getData(const char *auth_token, const char *dataset, float time, const char *spatial_interpolation, 
                 const char *temporal_interpolation, const char *variable, const char *spatial_operator, 
                 float points[][3], int num_points);
void writeVTK(const char *filename, double **data, int nx, int ny, int nz, int num_columns, float points[][3]);

#ifdef __cplusplus
}
#endif

#endif /* GETDATA_H */
