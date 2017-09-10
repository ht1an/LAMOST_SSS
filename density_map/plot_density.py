# plot the density map
# by yujc.astro@gmail.com
import sys
sys.path.append('../util/')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from aitoff_projection import aitoff_projection

dens = np.loadtxt('dens.dat')

## plot the outmost circle
#x0 = 0.0
#y0 = np.linspace(-90.0, 90.0, 18001)
#x1 = 360.0
#y1 = np.linspace(-90.0, 90.0, 18001)
#xm0, ym0 = np.meshgrid(np.deg2rad(x0), np.deg2rad(y0))
#xm1, ym1 = np.meshgrid(np.deg2rad(x1), np.deg2rad(y1))
#px0, py0 = aitoff_projection(xm0.T[0], ym0.T[0])
#px1, py1 = aitoff_projection(xm1.T[0], ym1.T[0])
##plt.plot(px0, py0, 'k-', px1, py1, 'k-', linewidth=1)


m = 3
nl = 360 * m
nb = 180 * m
delta = 1.0 / m
#l = np.linspace(  0.0, 360.0 - delta, nl)
#b = np.linspace(-90.0,  90.0 - delta, nb)
#l = np.linspace(  0.0 + delta/2, 360.0 - delta/2, nl)
#b = np.linspace(-90.0 + delta/2,  90.0 - delta/2, nb)
#d = np.zeros(nl*nb).reshape(nl, nb)
l = np.linspace(  0.0, 360.0, nl+1)
b = np.linspace(-90.0,  90.0, nb+1)
d = np.zeros((nl+1) * (nb+1)).reshape(nl+1, nb+1)
for i in range(360):
    for j in range(180):
        d[i*m:(i+1)*m, j*m:(j+1)*m] = dens[i, j]
d[nl, 0:nb] = d[nl-1, 0:nb]
d[0:nl, nb] = d[0:nl, nb-1]
d[nl, nb] = d[nl-1, nb-1]
lm, bm = np.meshgrid(np.deg2rad(l), np.deg2rad(b))
lm = lm.T
bm = bm.T
pl, pb = aitoff_projection(lm, bm) # projected l, b mesh
idx = (d > 0.0)
vmin = np.log10(np.min(d[idx]))
vmax = np.log10(np.max(d))

plt.pcolormesh(pl, pb, np.log10(d+1), vmin=vmin, vmax=vmax, cmap=cm.jet)
cbar = plt.colorbar()
plt.xlabel('l', fontdict={'fontsize':16})
plt.ylabel('b', fontdict={'fontsize':16})
cbar.set_label(r'$\log_{10}(\nu)$', fontdict={'fontsize':12})
#plt.savefig('gaia_densmap.png',format='png')
plt.show()
exit()
