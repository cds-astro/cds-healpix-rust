
# coding: utf-8

# In[28]:


import cdshealpix

from cdshealpix import cone_search_lonlat
from astropy_healpix.core import healpix_cone_search


# In[29]:


import mocpy

print("mocpy version: ", mocpy.__version__)
print("cdshealpix version: ", cdshealpix.__version__)


# In[33]:


import numpy as np
import astropy.units as u

lon = 20 * u.deg
lat = 0 * u.deg
radius = 50 * u.deg
depth = 2

# Get the cells from cdshealpix
cells_cds = cone_search_lonlat(lon=lon, lat=lat, radius=radius, depth=depth, delta_depth=1, flat=True)

# Get the cells from astropy_healpix
ipix_astropy = healpix_cone_search(lon=lon, lat=lat, radius=radius, nside=1 << depth, order="nested")
num_cells_astropy = ipix_astropy.shape[0]
cells_astropy = np.zeros(num_cells_astropy, dtype={
    'names':('ipix', 'depth', 'fully_covered'),
    'formats':(np.uint64, np.uint32, np.uint8),
})
cells_astropy["ipix"] = ipix_astropy
cells_astropy["depth"] = np.ones(num_cells_astropy) * depth
cells_astropy["fully_covered"] = np.ones(num_cells_astropy)


# In[34]:


from mocpy import MOC

# MOC instanciations
moc_cds = MOC.from_cells(cells_cds)
moc_astropy = MOC.from_cells(cells_astropy)


# In[36]:


import matplotlib.pyplot as plt
from astropy.wcs.utils import skycoord_to_pixel
from mocpy.spatial.utils import make_wcs

wcs = make_wcs(crpix=[0, 0], crval=[0, 0], cdelt=[-5, 5], ctype=["RA---AIT", "DEC--AIT"])

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=wcs)

# cdshealpix is green
moc_cds.fill(ax=ax,
             wcs=wcs,
             # mpl style kwargs
             edgecolor='g', facecolor='g', linewidth=1.0, fill=True, alpha=0.5)
# astropy_healpix is red
moc_astropy.fill(ax=ax,
                 wcs=wcs,
                 edgecolor='r', facecolor='r', linewidth=1.0, fill=True, alpha=0.25)

ax.axis('equal')
plt.xlabel('ra')
plt.ylabel('dec')
plt.grid(color='black', ls='dotted')
plt.title('cdshealpix: green, astropy_healpix: red')
plt.show()
plt.close()

# Check whether the cells returned by astropy_healpix are all contained in the cells returned by astropy_healpix
assert(np.in1d(cells_astropy["ipix"], cells_cds["ipix"]).all())

