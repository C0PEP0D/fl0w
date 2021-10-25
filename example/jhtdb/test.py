#!/usr/bin/python3

import os
import numpy as np
import pyJHTDB
import time

lJHTDB = pyJHTDB.libJHTDB()
lJHTDB.initialize()

data_set = 'isotropic1024coarse'
auth_token  = "com.gmail.remi.monthiller-fe721580"#"edu.jhu.pha.turbulence.testing-201311"  #Replace with your own token here
lJHTDB.add_token(auth_token)

print("Generating positions.")
time = 1.0
positions = np.zeros((1, 3), np.float32)
positions[0, :] = 1.0

print("Query positions.")
u = lJHTDB.getData(time, positions, sinterp=4, tinterp=0, data_set=data_set, getFunction='getVelocity')
j = lJHTDB.getData(time, positions, sinterp=40, tinterp=0, data_set=data_set, getFunction='getVelocityGradient')

print("flow.getVelocity(" + str(positions) + ", " + str(time) + ") -> " + str(u))
print("flow.getJacobian(" + str(positions) + ", " + str(time) + ") -> " + str(j.reshape((3,3))))

lJHTDB.finalize()

print("Start check.")
