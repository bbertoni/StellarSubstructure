from numpy import *
from astropy.io import fits
import re

#############################################################################

# import BOSS data:
hdulist = fits.open("/home/bridget/boss_data/bossstars5_cmr.fit")

objid = hdulist[1].data["objid"] # id number for the object
chunk = hdulist[1].data["chunk"] # survey chunk name (e.g. 'boss1')
lat = hdulist[1].data["glat"] # galactic latitude in deg
longi = hdulist[1].data["glong"] # galactic longitude in deg
z = hdulist[1].data["z"] # redshift to get LOS velocity
rmag = hdulist[1].data["psfMag_r"] # magnitude in the r band
imag = hdulist[1].data["psfMag_i"] # magnitude in the i band
gmag = hdulist[1].data["psfMag_g"] # magnitude in the g band
ext_r = hdulist[1].data["extinction_r"] # reddening in r band
ext_i = hdulist[1].data["extinction_i"] # reddening in i band
ext_g = hdulist[1].data["extinction_g"] # reddening in g band

# remove duplicate entries and only keep stars with abs(latitude) > 30 deg
icore, = where( (hdulist[1].data['sciencePrimary']==1) & (abs(lat) > 30) )

objidprim = objid[icore]
chunkprim = chunk[icore]
latprim = lat[icore]
longprim = longi[icore]
zprim = z[icore]
vprim = zprim * 2.99792458 * 10.**(5.0) # velocities in km/s
rmagprim = rmag[icore]
imagprim = imag[icore]
gmagprim = gmag[icore]
ext_rprim = ext_r[icore]
ext_iprim = ext_i[icore]
ext_gprim = ext_g[icore]

# get the object ids from the XDQSO data
objid_xdqso = []
gXDQSO = []
rXDQSO = []
iXDQSO = []
with open("/home/bridget/boss_data/XDQSO/mag_data/XDQSO_data_tot.txt",'r') as infile:
    for line in infile:
        params = line.split()
        objid_xdqso.append(int(params[0]))
        gXDQSO.append(float(params[1]))
        rXDQSO.append(float(params[2]))
        iXDQSO.append(float(params[3]))
objid_xdqso = array(objid_xdqso)
gXDQSO = array(gXDQSO)
rXDQSO = array(rXDQSO)
iXDQSO = array(iXDQSO)

# only pick out the stars with object id matches in the XDQSO dataset
idx = []
gXDQSONstars = []
rXDQSONstars = []
iXDQSONstars = []
for i in range(len(objid_xdqso)):
    g = where(objidprim==objid_xdqso[i])
    if size(g) != 0:
        idx.append(g[0][0])
        gXDQSONstars.append(gXDQSO[i])
        rXDQSONstars.append(rXDQSO[i])
        iXDQSONstars.append(iXDQSO[i])
gXDQSONstars = array(gXDQSONstars)
rXDQSONstars = array(rXDQSONstars)
iXDQSONstars = array(iXDQSONstars)

# final values for the XDQSO data set:

objidNstars = objidprim[idx]
chunkNstars = chunkprim[idx]
latNstars = latprim[idx] * pi/180.  # in radians
longNstars = longprim[idx] * pi/180.  # in radians
vNstars = vprim[idx]
rmagNstars = rmagprim[idx]
imagNstars = imagprim[idx]
gmagNstars = gmagprim[idx]
ext_rNstars = ext_rprim[idx]
ext_iNstars = ext_iprim[idx]
ext_gNstars = ext_gprim[idx]

N = len(latNstars)

cnumNstars = zeros(len(chunkNstars))

for i in range(len(chunkNstars)):
    imatch = re.search('[0-9]+',chunkNstars[i]) # matches any number in chunkNstars[i]
    cnumNstars[i] = int(imatch.group()) # retrieves the matches

##############################################################################

# find chunks for each star observation
cnum = zeros(len(chunk))

for i in range(len(chunk)):
    imatch = re.search('[0-9]+',chunk[i]) # matches any number in chunk[i]
    cnum[i] = int(imatch.group()) # retrieves the matches

# remove duplicate entries and only keep the stars in chunks that give rise to a uniform sample
# also only keep stars with abs(latitude) > 30 deg
icore2, = where((hdulist[1].data['sciencePrimary']==1) & (abs(lat) > 30) &
		(((cnum >= 14) &
		(hdulist[1].data['boss_target1'] & 2**40 != 0)) | (((cnum == 12) | (cnum == 13))
		& (hdulist[1].data['boss_target1'] & 2**40 != 0) & (hdulist[1].data['boss_target1'] 
		& 2**42 != 0))))
# 2**40 is QSO_CORE_MAIN, 2**42 is QSO_CORE_ED

# final values for the CORE data set:

objidCORE = objid[icore2]
chunkCORE = chunk[icore2]
latCORE = lat[icore2] * pi/180.  # in radians
longCORE = longi[icore2] * pi/180.  # in radians
zCORE = z[icore2]
vCORE = zCORE * 2.99792458 * 10.**(5.0) # velocities in km/s
rmagCORE = rmag[icore2]
imagCORE = imag[icore2]
gmagCORE = gmag[icore2]

hdulist.close()

##############################################################################

# convert velocities to the galactocentric rest frame
# coordinate system: x direction points from the sun to the galactic center,
# z direction points north, and right handed coordinates

# from arXiv:0912.3693 (velocity of sun wrt LSR):
u = 11.1 # positive x direction
v = 12.24 # positive y direction
w = 7.5 # positive z direction

vsun = 220.0 # positive y direction

vrf = ( vNstars + vsun*sin(longNstars)*cos(latNstars)
       + u*cos(longNstars)*cos(latNstars)
       + v*sin(longNstars)*cos(latNstars) + w*sin(latNstars))

################################################################################

# calculate the absolute magnitude from the r and i band colors using the bright
# normalization from http://iopscience.iop.org/article/10.1086/523619/pdf eqn. 1

Mr = ( 4.0 + 11.86*(rmagNstars - imagNstars)
       - 10.74*(rmagNstars - imagNstars)**2.0
       + 5.99*(rmagNstars - imagNstars)**3.0
       - 1.20*(rmagNstars - imagNstars)**4.0 )

# eqn. 2:
#Mr = ( 3.2 + 13.30*(rmagNstars - imagNstars)
#       - 11.50*(rmagNstars - imagNstars)**2.0
#       + 5.40*(rmagNstars - imagNstars)**3.0
#       - 0.70*(rmagNstars - imagNstars)**4.0 )

#print (N, len(Mr))
    
for i in range(N):
    if Mr[i] < 0:
        Mr[i] = NAN
    else:
        pass

################################################################################

# calculate the distance using eqn. 9 from
# http://iopscience.iop.org/article/10.1086/523619/pdf

dists = (10.0**( ((rmagNstars-Mr)/5.0 )+ 1.0)) * 0.001 # distance in kpc

################################################################################

# calculate directions and positions to each star

directions = zeros((N,3))
pos = zeros((N,3))

for i in range(N):
    directions[i,0] = cos(latNstars[i])*cos(longNstars[i])
    directions[i,1] = cos(latNstars[i])*sin(longNstars[i])
    directions[i,2] = sin(latNstars[i])
    pos[i,0] = dists[i]*directions[i,0]
    pos[i,1] = dists[i]*directions[i,1]
    pos[i,2] = dists[i]*directions[i,2]


###############################################################################

# output data file with positions, distance, LOS velocity (kpc and km/s),
# and chunk number

outfile = open("XDQSO_full_data.txt",'w')

for i in range(len(dists)):
    outfile.write(str(objidNstars[i])+" "+str(latNstars[i])+" "+str(longNstars[i])+
                  " "+str(vNstars[i])+" "+str(rmagNstars[i])+" "+str(imagNstars[i])+
                  " "+str(gmagNstars[i])+" "+str(ext_rNstars[i])+" "+str(ext_iNstars[i])+
                  " "+str(ext_gNstars[i])+" "+str(cnumNstars[i])+"\n")

outfile.close()

