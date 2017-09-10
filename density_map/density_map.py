# Calculate the number density per square degree for a given fits file
# by yujc.astro@gmail.com
import sys
import numpy as np
from astropy.io import fits

# our base catalogue, with the restriction that (dec >= -10, Gmag <= 15 deg)
# the table has only one extension
#input_fits = './base/base.fits'
n_argv = len(sys.argv)
if n_argv < 2:
    print('specific input fits name')
    sys.exit()
elif n_argv > 2:
    print('warning: more than one input fits provided, take the first one')
input_fits = str(sys.argv[1])
print('use file: ' + input_fits)

# try, throw, for furture update
hdulist = fits.open(input_fits)
tab = hdulist[1].data

l = tab['l']
b = tab['b']
ra = tab['ra']
dec = tab['dec']

# calculate the density use either ra-dec or galactic l-b coordinate
x = l
y = b
ymin = np.min(y)
ymax = np.max(y)
ymin = np.floor(ymin)
ymax = np.ceil(ymax)
deg2rad = np.pi / 180.0

# number density, density[i,j] refers to the density at (x = i, y = j - 90)
#density = np.zeros(360*180).reshape(360,180)
# number density, density[i,j] refers to the density at (x = j, y = i - 90)
density = np.zeros(360*180).reshape(180,360)

# number density calculated in an area larger than area_cr, 0.25 deg^2
area_cr = 0.25 
area_cr = area_cr * deg2rad # in unit of [deg * rad], facilitate the calculation
yi = int(ymax)
area_max = 0.0
while yi > ymin:
    dyi = 1
    dx = 1
    dy = np.sin(yi * deg2rad) - np.sin((yi-dyi) * deg2rad)
    if (dy >= area_cr):
        area = dx * dy
    else:
        if (dy * 360.0 < area_cr):
            dx = 360
            area = 0.0
            while (area < area_cr):
                dyi = dyi + 1
                dy = np.sin(yi * deg2rad) - np.sin((yi-dyi) * deg2rad)
                area = dx * dy
        elif (dy * 180.0 < area_cr):
            dx = 180
            area = 0.0
            while (area < area_cr):
                dyi = dyi + 1
                dy = np.sin(yi * deg2rad) - np.sin((yi-dyi) * deg2rad)
                area = dx * dy
        elif (dy * 90.0 < area_cr):
            dx = 90
            area = 0.0
            while (area < area_cr):
                dyi = dyi + 1
                dy = np.sin(yi * deg2rad) - np.sin((yi-dyi) * deg2rad)
                area = dx * dy
        elif (dy * 45.0 < area_cr):
            dx = 45
            area = 0.0
            while (area < area_cr):
                dyi = dyi + 1
                dy = np.sin(yi * deg2rad) - np.sin((yi-dyi) * deg2rad)
                area = dx * dy
        else:
            # max dx value: 5
            dx = 1
            area = 0.0
            while (area < area_cr):
                dyi = dyi + 1
                dy = np.sin(yi * deg2rad) - np.sin((yi-dyi) * deg2rad)
                dx = dx + 1
                if (dx > 5):
                    dx = 5
                area = dx * dy

    area = area / deg2rad
    if (yi - dyi < ymin):
        dyi = yi - ymin
        dy = np.sin(yi * deg2rad) - np.sin((yi-dyi) * deg2rad)
        area = dx * dy
    if (area > area_max):
        area_max = area

    # possible dx value: [1,2,3,4,5,45,90,180,360], so that 360%dx == 0
    idx = (y > yi - dyi) * (y <= yi)
    sub = x[idx]
    np.sort(sub)
    xi = 0
    dxi = dx
    while xi < 360:
        idx = (sub >= xi) * (sub < xi + dx)
        count = np.sum(idx)
        #density[xi:xi+dxi, 90+yi-dyi:90+yi] = count / area
        #print(xi,dxi,yi,dyi,count/area,density[xi,90+yi-dyi])
        density[90+yi-dyi:90+yi, xi:xi+dxi] = count / area
        #print(xi,dxi,yi,dyi,count/area,density[90+yi-dyi,xi])
        xi = xi + dxi
    yi = yi - dyi

np.savetxt('dens.dat',density.T)
