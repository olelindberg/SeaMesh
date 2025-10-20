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


# Remove the lines from the geo data frame:
gdf.drop(ids[1:],inplace=True)

#------------------------------------------------------------------------------
# Save:
#------------------------------------------------------------------------------
fp = project_root / "data/seamesh/Kystlinie_fixed.shp"
gdf.to_file(fp)

#------------------------------------------------------------------------------
# Plot:
#------------------------------------------------------------------------------
for geom in gdf.geometry:
    geom_type = geom.geom_type
    print(f"Geometry type: {geom_type}")
    x,y = geom.xy
    plt.plot(x,y)
plt.show()
