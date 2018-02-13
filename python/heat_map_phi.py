from __future__ import print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap as Basemap
from matplotlib.colors import rgb2hex
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import csv

df = pd.read_csv('4GRP_spatial_difference_summary_rev.csv', usecols=[
                 'state', 'abb', 'rho', 'phi', 'theta', 'total', 'rho_per', 'phi_per', 'theta_per', 'latitude', 'longitude'])

df_location = pd.read_csv('location.csv', usecols=['state', 'abb', 'capital', 'latitude', 'longitude'])

# Lambert Conformal map of lower 48 states.
plt.figure()
m = Basemap(llcrnrlon=-119,llcrnrlat=22,urcrnrlon=-64,urcrnrlat=49,
            projection='lcc',lat_1=33,lat_2=45,lon_0=-95)

shp_info = m.readshapefile('st99_d00','states',drawbounds=True)
# shp_info = m.readshapefile('gz_2010_us_040_00_500k', 'states',drawbounds=True)

# Dictionary {'yyy':xxx, ..., }
# dic.keys() -> Return the index
drho ={k: float(v) for k,v in df.groupby("state")["rho"]}
dphi ={k: float(v) for k,v in df.groupby("state")["phi"]}
dtheta ={k: float(v) for k,v in df.groupby("state")["theta"]}

A = df['phi'].values

exclusions = ['District of Columbia','Puerto Rico', 'Rhode Island', 'Delaware', 'Hawaii', 'New Hampshire', 'South Carolina','Louisiana', 'Kentucky', 'Missouri', 'West Virginia', 'Vermont', 'Oklahoma', 'Arizona', 'Maine', 'Kansas', 'Utah', 'Nebraska', 'Nevada', 'Idaho', 'New Mexico', 'South Dakota', 'North Dakota', 'Montana','Wyoming','Alaska']

lons = df['longitude'].values
lats = df['latitude'].values
abbs = df['abb'].values
percentages = df['phi_per'].values
include = df['state'].values
# lons =[]
# lats = []
# for i, statename in enumerate(df['state']):
#     lons.append(df['longitude'][i])
#     lats.append(df['latitude'][i])


# for i, statename in enumerate(df_location['state']):
#     if statename not in exclusions: 
#         lons.append(df_location['longitude'][i])
#         lats.append(df_location['latitude'][i])

print(shp_info)
# choose a color for each state based on population density.
colors={}
statenames=[]
patches = []
color_values =[]

cmap = plt.cm.jet # use 'hot' colormap plt.cm.bwr #
vmin = 1.0; vmax = 2.0 # set range.
phi_max = max(A)
phi_min = min(A)
print(m.states_info[0].keys())

for shapedict in m.states_info:
    statename = shapedict['NAME']
    # skip DC and Puerto Rico.
    if statename in include: #not in exclusions: #['District of Columbia','Puerto Rico']:
        phis = dphi[statename]
        # calling colormap with value between 0 and 1 returns
        # rgba value.  Invert color range (hot colors are high
        # population), take sqrt root to spread out colors more.
        # color_value = 1.-np.sqrt((phis-vmin)/(vmax-vmin))
        color_value = (phis-phi_min)/(phi_max-phi_min)
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
p = PatchCollection(patches, cmap=plt.cm.jet, alpha=0.6)
# p.set_array(np.array(A))
# p.set_array(np.array(color_values))
re_color = np.array(color_values) * (phi_max-phi_min)  + phi_min
p.set_array(re_color)
# p.set_array(np.array(colors))
ax.add_collection(p)
# plt.colorbar(p, orientation="horizontal")
plt.colorbar(p, shrink=0.7,orientation="horizontal")


# c = plt.colorbar(cmap, orientation='horizontal')
# fig  = plt.gci()
# m.colorbar()
# cb = m.colorbar("bottom", size="5%", pad="2%")
# draw meridians and parallels.
# m.drawparallels(np.arange(25,65,20),labels=[1,0,0,0])
# m.drawmeridians(np.arange(-120,-40,20),labels=[0,0,0,1])
# plt.title('Filling State Polygons by Population Density')
for lon, lat, abb, percentage  in zip(lons, lats, abbs, percentages):
    x,y = m(lon, lat)
    # msize = mag * min_marker_size
    # marker_string = get_marker_color(mag)
    # m.plot(x, y, 'bo', markersize=6)
    # state = df_location['state'] if df_location['latitude'] == lat
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
        position = (10,-20)
    
    if abb not in ['MA', 'CT', 'NJ', 'MD']:
        plt.annotate('%s \n %d %%' % (abb, percentage) , 
            xy = (x, y), xytext = position, 
            textcoords = 'offset points', ha = 'right', va = 'bottom')
    else:
        m.plot(x, y, 'ro', markersize=6)
        plt.annotate('%s \n %d %%' % (abb, percentage) , 
            xy = (x, y), xytext = position, 
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'None', alpha = 0.5),
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

plt.show()
