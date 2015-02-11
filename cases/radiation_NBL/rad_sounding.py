import numpy as np
from pylab import *

# ------------------
# Settings
# ------------------
ps   = 1e5     # surface pressure
th0  = 290     # initial mixed-layer temperature
RH0  = 0.5     # initial (lowest model level) relative humidity
ztop = 1000    # estimated domain top
dz0  = 10       # vertical grid spacing lowest level

# ------------------
# make stretched grid
# ------------------
# 1. Calculate at half levels for zm_grid_in
zh = ([0])
dz = ([dz0])
while(zh[-1] < ztop):
    if(zh[-1] < 200):
        dzrat = 1.
    elif(zh[-1] < 400):
        dzrat = 1.03
    else:
        dzrat = 1.1
            
    dz0 *= dzrat
    zh.append(zh[-1]+dz0)
    dz.append(dz0)

# Cast to numpy array
zh = np.array(zh)
dz = np.array(dz)

# 2. Calculate full levels:
z = 0.5 * (zh[1:] + zh[:-1])

print('%i half levels in zm_grid_in'%zh.size)

# -------------------
# make sounding
# -------------------
thl = np.zeros(z.size)
qt  = np.zeros(z.size)
u   = np.zeros(z.size)
v   = np.zeros(z.size)

thl[:] = th0
u  [:] = 5.
v  [:] = 0.

# Calculate moisture mixing ratio
es = 0.611e3 * np.exp(17.2694 * (thl[0] - 273.16) / (thl[0] - 35.86))
qs = 0.622 * es / ps
qt[:] = RH0 * qs

# -------------------
# Write to file
# -------------------
gridfile = open('zm_grid_in', 'w')
for k in range(z.size):
    gridfile.write('%14.8e\n'%(zh[k]))
gridfile.close()

proffile = open('sound_in', 'w')
proffile.write('%14.8e %14.8e %14.8e %14.8e %14.8e\n'%(ps/100., thl[0], qt[0]*1000, u[0], v[0]))
for k in range(z.size):
    proffile.write('%14.8e %14.8e %14.8e %14.8e %14.8e\n'%(z[k], thl[k], qt[k]*1000, u[k], v[k]))
proffile.close()

# ------------------
# Plot, including backrad sounding
# ------------------
br = np.loadtxt('backrad_in', skiprows=1)
bp = br[:,0][::-1]
bt = br[:,1][::-1]
bq = br[:,2][::-1]



figure()
subplot(121)
plot(bt,bp)
subplot(122)
plot(bq,bp)
