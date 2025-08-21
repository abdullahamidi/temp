#!/usr/bin/env /bin/python3

import numpy as np
import math
from scipy.interpolate import Rbf
import sys
import os.path
scriptDir = os.path.dirname(__file__);

# Yates
massConv = 175.126     # 1 lbf*in/s^2 = 175.125 kg
lenConv  = 0.0254      # 1 inch = 0.0254 m

damping  = 0.02        # Experiment
#damping  = 0.0         # FOI (Edge) analysis assume no structural damping


## Mass, damping and stiffness matrices ................................

nMS      = 3           # Number of mode shapes
M        = [[0.0 for x in range(nMS)] for y in range(nMS)]
f        = [[0.0 for x in range(nMS)] for y in range(nMS)]
w        = [[0.0 for x in range(nMS)] for y in range(nMS)]
C        = [[0.0 for x in range(nMS)] for y in range(nMS)]
K        = [[0.0 for x in range(nMS)] for y in range(nMS)]

df     = np.loadtxt(os.path.join(scriptDir, "yates/mass"), comments = '#')
mode   = df[:,0]
mass   = df[:,1]
df     = np.loadtxt(os.path.join(scriptDir, "yates/frequency"), comments = '#')
freq   = df[:,1]

for i in range(nMS):
    M[i][i] = massConv*mass[i]
    f[i][i] = freq[i]

w      = np.dot(2.0*np.pi,f)            # Angular frequency - 2 pi f
K      = np.dot(w,np.dot(w,M))          # Stiffness         - w w M
C      = np.dot(2*damping,np.dot(w,M))  # Damping           - 2 xi w M


# Read mode shape .....................................................
df     = np.loadtxt("points", comments = '#')
x      = df[:,0]
y      = df[:,1]
z      = df[:,2]


# Write mass, damping and stiffness matrices ..............................

orig_stdout = sys.stdout

fOut = 'massMatrix' 
f = open(fOut, 'w')
sys.stdout = f

print("massMatrix "+str(nMS*nMS)+"(")
for i in range(nMS):
    for j in range(nMS):
        sys.stdout.write(str(M[i][j])+" ")
    print("")
print(");")

sys.stdout = orig_stdout
f.close()


orig_stdout = sys.stdout

fOut = 'dampingMatrix'
f = open(fOut, 'w')
sys.stdout = f

print("dampingMatrix "+str(nMS*nMS)+"(")
for i in range(nMS):
    for j in range(nMS):
        sys.stdout.write(str(C[i][j])+" ")
    print("")
print(");")

sys.stdout = orig_stdout
f.close()


orig_stdout = sys.stdout

fOut = 'stiffnessMatrix'
f = open(fOut, 'w')
sys.stdout = f

print("stiffnessMatrix "+str(nMS*nMS)+"(")
for i in range(nMS):
    for j in range(nMS):
        sys.stdout.write(str(K[i][j])+" ")
    print("")
print(");")

sys.stdout = orig_stdout
f.close()


# Yates mode shapes .......................................................

df    = np.loadtxt(os.path.join(scriptDir, "yates/vertices"), comments = '#')
xE    = lenConv*df[:,1]
yE    = lenConv*df[:,2]
zE    = lenConv*df[:,3]

# Mode shapes are dimensionless
df    = np.loadtxt(os.path.join(scriptDir, "yates/mode1"), comments = '#')
xM1E  = df[:,1]
yM1E  = df[:,2]
zM1E  = df[:,3]

df    = np.loadtxt(os.path.join(scriptDir, "yates/mode2"), comments = '#')
xM2E  = df[:,1]
yM2E  = df[:,2]
zM2E  = df[:,3]

df    = np.loadtxt(os.path.join(scriptDir, "yates/mode3"), comments = '#')
xM3E  = df[:,1]
yM3E  = df[:,2]
zM3E  = df[:,3]


# Interpolate measured values to patch points .............................

rbfM1x = Rbf(xE, yE, xM1E, function='linear')
rbfM1y = Rbf(xE, yE, yM1E, function='linear')
rbfM1z = Rbf(xE, yE, zM1E, function='linear')
rbfM2x = Rbf(xE, yE, xM2E, function='linear')
rbfM2y = Rbf(xE, yE, yM2E, function='linear')
rbfM2z = Rbf(xE, yE, zM2E, function='linear')
rbfM3x = Rbf(xE, yE, xM3E, function='linear')
rbfM3y = Rbf(xE, yE, yM3E, function='linear')
rbfM3z = Rbf(xE, yE, zM3E, function='linear')

xM1    = rbfM1x(x,y)
yM1    = rbfM1y(x,y)
zM1    = rbfM1z(x,y)
xM2    = rbfM2x(x,y)
yM2    = rbfM2y(x,y)
zM2    = rbfM2z(x,y)
xM3    = rbfM3x(x,y)
yM3    = rbfM3y(x,y)
zM3    = rbfM3z(x,y)



## Write mode shape ........................................................

orig_stdout = sys.stdout

fOut = 'modeShapes' 
f = open(fOut, 'w')
sys.stdout = f

print("modeShapes")
print(int(nMS + (nMS+1.0)*nMS/2.0))
print("( \n"+str(len(x))+"\n(")

# Mode 1
for i in range(len(x)):
    print("(%.8e %.8e %.8e)" % (xM1[i], yM1[i], zM1[i]))

print(") \n\n"+str(len(x))+"\n(")

# Mode 2
for i in range(len(x)):
    print("(%.8e %.8e %.8e)" % (xM2[i], yM2[i], zM2[i]))

print(") \n\n"+str(len(x))+"\n(")

# Mode 3
for i in range(len(x)):
    print("(%.8e %.8e %.8e)" % (xM3[i], yM3[i], zM3[i]))

print(")")

# Quadratic modes are neglected for this analysis
print(" ")
print(str(len(x))+"{(0.0 0.0 0.0)}\n")   # sigma_11
print(str(len(x))+"{(0.0 0.0 0.0)}\n")   # sigma_12
print(str(len(x))+"{(0.0 0.0 0.0)}\n")   # sigma_13
print(str(len(x))+"{(0.0 0.0 0.0)}\n")   # sigma_22
print(str(len(x))+"{(0.0 0.0 0.0)}\n")   # sigma_23
print(str(len(x))+"{(0.0 0.0 0.0)}\n")   # sigma_33

print(");")

sys.stdout = orig_stdout
f.close()
