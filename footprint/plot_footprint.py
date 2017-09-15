# plot the density map
# by yujc.astro@gmail.com
import sys
sys.path.append('../util/')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from aitoff_projection import aitoff_projection
from footprint import footprint
import galactic

fig, ax = plt.subplots(1,1)
# plot coordinate
# x, y = ra, dec
dx = 60
x  = 0
y0 = np.linspace(-90.0, 90.0, 18001)
while x <= 360:
    x0 = np.array([x]*len(y0))
    px0, py0 = aitoff_projection(np.deg2rad(x0), np.deg2rad(y0))
    if x == 0 or x == 360:
        ax.plot(px0, py0, 'k-', lw=2)
    elif x == 180:
        ax.plot(px0, py0, 'k-', lw=1)
    else:
        ax.plot(px0, py0, 'k--', lw=1)
    x = x + dx

x0 = np.linspace(0.0, 360.0, 36001)
dy = 30
y = -90 + dy
while y < 90:
    y0 = np.array([y]*len(x0))
    px0, py0 = aitoff_projection(np.deg2rad(x0), np.deg2rad(y0))
    if y == 0:
        ax.plot(px0, py0, 'k-', linewidth=1)
    else:
        ax.plot(px0, py0, 'k--', linewidth=1)
    y = y + dy

# plot galactic latitude = 30, 60
b = 0
db = 30
l0 = np.linspace(0.0, 360.0, 36001)
while b < 90:
    b0 = np.array([b]*len(l0))
#    c = SkyCoord(l0, b0, unit='deg', frame='galactic')
#    x0 = c.icrs.ra.value
#    y0 = c.icrs.dec.value

    x0, y0 = galactic.galactic2(l0, b0)
    if b == 0:
        ra, dec = galactic.galactic2(l0, b0)
        idx = np.argsort(ra)
        x0 = ra[idx]
        y0 = dec[idx]
    else:
        x0, y0 = galactic.galactic2(l0, b0)
    px0, py0 = aitoff_projection(np.deg2rad(x0), np.deg2rad(y0))
    ax.plot(px0, py0, 'b-.', linewidth=0.5, alpha=0.5)
    b = b + db

gp_ra = 192.859508
gp_dec = 27.128336
px1, py1 = aitoff_projection(np.deg2rad(gp_ra), np.deg2rad(gp_dec))
ax.scatter(px1, py1, s=4, c='b', marker='*', alpha=0.5)

ax.text(-12,  -3, r'0$^{\circ}$',   fontdict={'fontsize':8})
ax.text( 11,  45, r'30$^{\circ}$',  fontdict={'fontsize':8})
ax.text( 76,  78, r'60$^{\circ}$',  fontdict={'fontsize':8})
ax.text( 36,  -6, r'60$^{\circ}$',  fontdict={'fontsize':8})
ax.text( 92,  -6, r'120$^{\circ}$', fontdict={'fontsize':8})
ax.text(158,  -6, r'180$^{\circ}$', fontdict={'fontsize':8})
ax.text(224,  -6, r'240$^{\circ}$', fontdict={'fontsize':8})
ax.text(285,  -6, r'300$^{\circ}$', fontdict={'fontsize':8})
ax.text( 7,  -51, r'-30$^{\circ}$', fontdict={'fontsize':8})
ax.text( 72, -83, r'-60$^{\circ}$', fontdict={'fontsize':8})
ax.text(195,  25, r'NGP', fontdict={'fontsize':8}, color='red', alpha=0.3)

t = np.genfromtxt('../script/center.txt',dtype=('S19', float, float, 'S19'), names=True)
c1 = t['ra']
c2 = t['dec']
name1 = t['name']
name2 = np.array(list(set(name1)))

colors = cm.rainbow(np.linspace(0, 1, len(name2)))

#for nx in name2:
for j in range(len(name2)):
    nx = name2[j]
    idx = (name1 == nx)
    cen_ra  = c1[idx]
    cen_dec = c2[idx]
    if (nx == 'JiangDengKai'):
        k = j
        continue
    for i in range(len(cen_ra)):
        x0, y0 = footprint(cen_ra[i], cen_dec[i])
        x, y = aitoff_projection(np.deg2rad(x0), np.deg2rad(y0))
        if i == 0:
            ax.plot(x, y, c=colors[j], ls='-', label=nx)
        else:
            ax.plot(x, y, c=colors[j], ls='-')

nx = name2[k]
idx = (name1 == nx)
cen_ra  = c1[idx]
cen_dec = c2[idx]
x0, y0 = footprint(cen_ra, cen_dec)
x, y = aitoff_projection(np.deg2rad(x0), np.deg2rad(y0))
ax.plot(x, y, c=colors[k], ls='-', label=nx)
ax.legend(loc=1, bbox_to_anchor=(1.06, 1.06), fontsize=8, fancybox=True, framealpha=0.5)

ax.set_title('LAMOST-MRS test footprint')

plt.savefig('fp.pdf', format='pdf')
