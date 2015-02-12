import numpy as np
from pylab import *
from netCDF4 import Dataset

# ------------------
# Settings
# ------------------
dz0   = 25.   ; dzrat = 1.03
#dz0   = 12.5  ; dzrat = 1.06
#dz0   = 6.25  ; dzrat = 1.07
#dz0   = 3.125 ; dzrat = 1.075
#dz0   = 2.    ; dzrat = 1.08
#dz0   = 1.    ; dzrat = 1.08

ps    = 1e5     # surface pressure
th0   = 290     # initial mixed-layer temperature
zi0   = 900
dthdz = 7e-3    # temperature gradient free tropospheredqdz  = -2e-6   # moisture gradient free troposphere
dqdz  = -1.3e-6   # moisture gradient free troposphere
RH0   = 0.52    # initial (lowest model level) relative humidity
ztop  = 1500    # estimated domain top
maxdz = 40      # maximum vertical grid spacing

# ------------------
# make stretched grid
# ------------------
# 1. Calculate at half levels for zm_grid_in
zh = ([0])
dz = ([dz0])
while(zh[-1] < ztop):
    if(zh[-1] < 400):
        rat = 1.
    else:
        rat = dzrat
            
    dz0 *= rat
    dz0 = np.min((dz0, maxdz))

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

# Calculate moisture mixing ratio
es  = 0.611e3 * np.exp(17.2694 * (th0 - 273.16) / (th0 - 35.86))  # only for ps=pref
qs  = 0.622 * es / ps
qt0 = RH0 * qs

u  [:] = 7.
v  [:] = 0.
for k in range(z.size):
    if(z[k] < zi0):
        thl[k] = th0
        qt [k] = qt0
    else:
        thl[k] = th0 + (z[k]-zi0) * dthdz
        qt [k] = qt0 + (z[k]-zi0) * dqdz

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
if(False):
    #close('all')

    figure()
    subplot(221)
    plot(thl,z)
    subplot(222)
    plot(qt,z)
    subplot(223)
    plot(dz,zh)


if(True):
    close('all')

    # Read backrad_in
    br = np.loadtxt('backrad_in', skiprows=1)
    bp = br[:,0][::-1]
    bt = br[:,1][::-1]
    bq = br[:,2][::-1]
    bth = bt / (bp / 1e3)**(287/1004)
   
    # Read output UCLA to get correct pressure, density, etc.
    nc = Dataset('rad.ps.nc','r')
    nzt = nc.variables["zt"][:]
    np  = nc.variables["p"][0,:]
    nth = nc.variables["t"][0,:]
    nr  = nc.variables["q"][0,:]
    nT  = nth * (np / 1e5)**(287/1004)
    
    figure()
    subplot(121)
    plot(bth,bp*100,label='backrad')
    plot(nth,np)
    xlim(280,310)
    ylim(60000,101300)

    subplot(122)
    plot(bq*1000,bp*100,label='backrad')
    plot(nr,np)
    #xlim(280,310)
    ylim(60000,101300)
