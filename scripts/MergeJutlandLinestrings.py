# Import necessary modules
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import shapely as shp
from shapely.geometry import Polygon
from pathlib import Path


project_root = Path("/home/ole/Projects/SeaMesh/")
fp           = project_root / "data/Topo-bathy/shape_format/data/Kystlinie.shp"
gdf          = gpd.read_file(fp)

#------------------------------------------------------------------------------
# Fix: connect linestrings in Jylland:
#------------------------------------------------------------------------------
ids = []
for i1, geom1 in enumerate(gdf.geometry):
    for i2, geom2 in enumerate(gdf.geometry):
        if (geom1.touches(geom2)):

            if not i1 in ids:
                ids.append(i1)

            if not i2 in ids:
                ids.append(i2)

# Merge the lines together:
line = gdf.geometry[ids[0]]
for id in ids[1:]:
    line = shp.ops.linemerge([line,gdf.geometry[id]])

# Assign the new line:
gdf.at[ids[0],"geometry"] = line


#------------------------------------------------------------------------------
# Convert linestrings to polygons
#------------------------------------------------------------------------------
gdf["geometry"] = gdf["geometry"].apply(lambda geom: Polygon(geom))

#------------------------------------------------------------------------------
# Outer bounding box boundary:
#------------------------------------------------------------------------------
minx, miny, maxx, maxy = gdf.total_bounds
width                  = min(maxx-minx,maxy-miny)/2
bbox_polygon           = shp.geometry.box(minx-width, miny-width, maxx+width, maxy+width)
gdf.at[ids[1],"geometry"] = bbox_polygon

# Remove the lines from the geo data frame:
gdf.drop(ids[2:],inplace=True)

#------------------------------------------------------------------------------
# Save:
#------------------------------------------------------------------------------
fp = project_root / "data/seamesh/Kystlinie_fixed.shp"
gdf.to_file(fp)

#------------------------------------------------------------------------------
# Plot:
#------------------------------------------------------------------------------
for polygon in gdf.geometry:
    x,y = polygon.exterior.xy
    plt.plot(x,y)
plt.show()
