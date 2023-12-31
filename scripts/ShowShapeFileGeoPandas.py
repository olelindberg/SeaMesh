# Import necessary modules
import geopandas as gpd
import matplotlib.pyplot as plt

# Set filepath
fp = "../data/Topo-bathy/shape_format/data/Kystlinie.shp"
# fp = "../data/Topo-bathy/shape_format/data/Hav_dybde_5m_2005.shp"
# fp = "../data/Topo-bathy/shape_format/data/Land_hojde_5m_2005.shp"
# fp = "../data/coastlines-split-3857/lines.shp"
# fp = "../data/coastlines-split-4326/lines.shp"

# Read file using gpd.read_file()
data = gpd.read_file(fp)
type(data)

print(data["geometry"][0])
print(data.head(2))

data.plot()
plt.show()
