from __future__ import print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap as Basemap
from matplotlib.colors import rgb2hex
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import Normalize
import csv

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))


df = pd.read_csv('diff_social_planner.csv', usecols=[
                 'state', 'abb', 'nu','welfare','yk', 'ok', 'yl', 'ol', 'latitude', 'longitude'])

df_location = pd.read_csv('location.csv', usecols=['state', 'abb', 'capital', 'latitude', 'longitude'])

# Lambert Conformal map of lower 48 states.
plt.figure()
m = Basemap(llcrnrlon=-119,llcrnrlat=22,urcrnrlon=-64,urcrnrlat=49,
            projection='lcc',lat_1=33,lat_2=45,lon_0=-95)

shp_info = m.readshapefile('st99_d00','states',drawbounds=True)
# shp_info = m.readshapefile('gz_2010_us_040_00_500k', 'states',drawbounds=True)

# Dictionary {'yyy':xxx, ..., }
# dic.keys() -> Return the index
dX ={k: float(v) for k,v in df.groupby("state")["ol"]}
A = df['ol'].values

exclusions = ['District of Columbia','Puerto Rico', 'Rhode Island', 'Delaware', 'Hawaii', 'New Hampshire', 'South Carolina','Louisiana', 'Kentucky', 'Missouri', 'West Virginia', 'Vermont', 'Oklahoma', 'Arizona', 'Maine', 'Kansas', 'Utah', 'Nebraska', 'Nevada', 'Idaho', 'New Mexico', 'South Dakota', 'North Dakota', 'Montana','Wyoming','Alaska']

lons = df['longitude'].values
lats = df['latitude'].values
abbs = df['abb'].values
# percentages = df['phi_per'].values
include = df['state'].values

print(shp_info)
# choose a color for each state based on population density.
colors={}
statenames=[]
patches = []
color_values =[]

cmap = plt.cm.seismic #cool # use 'hot' colormap plt.cm.bwr #
# vmin = 0.06; vmax = 0.12 # set range.

X_max = 0.3
X_min = -0.1
print(m.states_info[0].keys())

for shapedict in m.states_info:
    statename = shapedict['NAME']
    # skip DC and Puerto Rico.
    if statename in include: #not in exclusions: #['District of Columbia','Puerto Rico']:
        Xs = dX[statename]
        # calling colormap with value between 0 and 1 returns
        # rgba value.  Invert color range (hot colors are high
        # population), take sqrt root to spread out colors more.
        # color_value = 1.-np.sqrt((Xs-vmin)/(vmax-vmin))
        # if statename == 'New York':
        #     Xs = X_max
        # if statename == 'Pennsylvania':
        #     Xs = X_max
        # if statename == 'Massachusetts':
        #     Xs = X_min
        color_value = (Xs-X_min)/(X_max-X_min)
        colors[statename] = cmap(color_value)[:3]
        color_values.append(color_value)
    statenames.append(statename)
# cycle through state names, color each one.

ax = plt.gca() # get current axes instance
polys = np.array([])

for nshape,seg in enumerate(m.states):
    # skip DC and Puerto Rico.
    if statenames[nshape] not in exclusions: #['District of Columbia','Puerto Rico']:
        color = rgb2hex(colors[statenames[nshape]]) 
        poly = Polygon(seg,facecolor=color,edgecolor=color)
        # ax.add_patch(poly)
        patches.append(poly)
  
# colors = np.random.rand(len(patches))
norm = MidpointNormalize(midpoint=0)
p = PatchCollection(patches, norm=norm, cmap=plt.cm.seismic, alpha=0.8)
# p.set_array(np.array(A))
# p.set_array(np.array(color_values))
re_color = np.array(color_values) * (X_max-X_min)  + X_min
p.set_array(re_color)
# p.set_array(np.array(colors))
ax.add_collection(p)
colorbar = plt.colorbar(p, shrink=0.9, orientation="horizontal")
mn = 0.3
mx = -0.1
md = 3
colorbar.set_ticks(np.linspace(mn, mx, md))
# colorbar.set_ticklabels(np.linspace(mn, mx, md))

f_size = 15

for lon, lat, abb  in zip(lons, lats, abbs):
    x,y = m(lon, lat)
    if abb=='WA':
        position = (35, -10)
    elif abb=='OR':
        position = (20, -25)
    elif abb=='CA':
        position = (30, -60)
    elif abb=='MN':
        position = (0, 0)
    elif abb=='IA':
        position = (10, -10)
    elif abb=='WI':
        position = (10, 5)
    elif abb=='MS':
        position = (15, -20)
    elif abb=='GA':
        position = (26, -25)
    elif abb=='FL':
        position = (10, -40)
    elif abb=='IL':
        position = (16, -10)
    elif abb=='IN':
        position = (14, -10)
    elif abb=='OH':
        position = (10, -10)
    elif abb=='MI':
        position = (12, -15)
    elif abb=='NY':
        position = (-2, -5)
    elif abb=='PA':
        position = (-5, -10)
    elif abb=='VA':
        position = (0, -15)
    elif abb=='MA':
        position = (-10, 50)
    elif abb=='CT':
        position = (30, -70)
    elif abb=='NJ':
        position = (45, -90)
    elif abb=='MD':
        position = (30, -120)
    else:
        position = (10,-10)
    
    if abb not in ['MA', 'CT', 'NJ', 'MD']:
        plt.annotate('%s' % abb , 
            xy = (x, y), xytext = position, fontsize=f_size,
            textcoords = 'offset points', ha = 'right', va = 'bottom')
    else:
        m.plot(x, y, 'ro', markersize=6)
        plt.annotate('%s' % abb, 
            xy = (x, y), xytext = position, fontsize=f_size,
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'None', alpha = 0.5),
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

plt.show()
