from numpy import *
from astropy.io import fits
import sys
import re
from matplotlib import pyplot as plt

#############################################################################

# import XDQSO data:

infile = open("XDQSO_vel_corr_data.txt",'r')

dists = []
vrf = []
cnum = []

for line in infile:
	params = line.split()
	dists.append(float(params[3]))
	vrf.append(float(params[4]))
	cnum.append(int(params[5]))

dists = array(dists)
vrf = array(vrf)
cnum = array(cnum)

dists = dists[cnum > 1]
vrf = vrf[cnum > 1]

dists = dists[dists > 5]
vrf = vrf[dists > 5]

###############################################################################

# plot results

fig = plt.figure()
plt.plot(dists[abs(vrf)<600],vrf[abs(vrf)<600],'b.',alpha=0.2,markersize=4)
plt.xlim(0,300)
plt.xlabel("distance from sun (kpc)",size=14)
plt.ylabel("line of sight velocity (km/s)",size=14)
plt.title("XDQSO data, |b| < 30 deg, vLOS < 600 km/s",size=15)

fig.set_size_inches(23.21, 12.42)

plt.show()

fig.savefig("XDQSO_vLOS_dist_dcut.png",dpi=100)

