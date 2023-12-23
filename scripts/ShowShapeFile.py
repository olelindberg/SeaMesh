import matplotlib.pyplot as plt
import shapefile as sf
import numpy as np

filepath = "../data/Topo-bathy/shape_format/data/Kystlinie.shp"
# filepath = "../data/Topo-bathy/shape_format/data/Hav_dybde_5m_2005.shp"
# filepath = "../data/Topo-bathy/shape_format/data/Land_hojde_5m_2005.shp"
# filepath = "../data/coastlines-split-3857/lines.shp"
# filepath = "../data/coastlines-split-4326/lines.shp"

# Read file using gpd.read_file()
shapefile = sf.Reader(filepath)

bbox = np.reshape(np.array(shapefile.bbox), (2, 2))
print(np.reshape(bbox, (2, 2)))
shapes = shapefile.shapes()


plt.plot(bbox[:, 0], [bbox[0, 1], bbox[0, 1]])
plt.plot(bbox[:, 0], [bbox[1, 1], bbox[1, 1]])
plt.plot([bbox[0, 0], bbox[0, 0]], bbox[:, 1])
plt.plot([bbox[1, 0], bbox[1, 0]], bbox[:, 1])

for shape in shapes:
    points = np.array(shape.points)
    plt.plot(points[:, 0], points[:, 1])


plt.show()
