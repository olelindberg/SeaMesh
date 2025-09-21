import matplotlib.pyplot    as plt
import numpy                as np
import triangle             as tr
import shapely              as shp
import geopandas            as gpd

from pathlib import Path

#------------------------------------------------------------------------------
# Load shapefile with polygons:
#------------------------------------------------------------------------------
project_root = Path("/home/ole/Projects/SeaMesh/")
fp           = project_root / "data/seamesh/Kystlinie_fixed.shp"
gdf          = gpd.read_file(fp)
gdf          = gdf.assign(area=gdf.area).sort_values("area",ascending=False)


#------------------------------------------------------------------------------
# Extract vertices and line segments from polygons:
#------------------------------------------------------------------------------
vertices_list       = []
linesegments_list   = []
points_inside       = []
offset              = 0
for poly_id, polygon in enumerate(gdf.geometry):
    
    # Find a point inside the polygon:
    candidate = polygon.centroid
    if polygon.contains(candidate):
        p = candidate
    else:
        # fallback: just pick first interior point from buffer trick
        p = polygon.representative_point()
    points_inside.append([p.x,p.y])

    # Get polygon vertices:
    x,y = polygon.exterior.xy

    # Add vertices (skip the last point since it is the same as the first):
    vertices = np.stack([x[:-1], y[:-1]], axis=1)
    vertices_list.append(vertices)

    # Add line segments:
    num_pts      = shp.get_num_coordinates(polygon)-1
    i            = np.arange(num_pts)
    linesegments = np.stack([i, i + 1], axis=1) % num_pts
    linesegments_list.append(linesegments+offset)
    
    # Increment offset:
    offset = offset + num_pts

#------------------------------------------------------------------------------
# Triangulate with triangle:
#------------------------------------------------------------------------------
vertices        = np.vstack(vertices_list)
linesegments    = np.vstack(linesegments_list)
tridata         = dict(vertices=vertices, segments=linesegments,holes=points_inside[1:])#, holes=[[0, 0]])
trimesh         = tr.triangulate(tridata,'q30p')
vertices        = trimesh["vertices"]
triangles       = trimesh["triangles"]
print(f"Triangle mesh with {vertices.shape[0]} vertices and {triangles.shape[0]} triangles")

#------------------------------------------------------------------------------
# Plot:
#------------------------------------------------------------------------------
plt.triplot(vertices[:,0],vertices[:,1],triangles)
for polygon in gdf.geometry:
    x,y = polygon.exterior.xy
    plt.plot(x,y)
plt.axis("equal")
plt.show()