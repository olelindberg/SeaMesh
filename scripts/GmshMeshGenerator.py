import geopandas as gpd
import gmsh 
import sys
from pathlib import Path
length_scale = 100000

project_root = Path("/home/ole/Projects/SeaMesh/")
fp           = project_root / "data/seamesh/Kystlinie_fixed.shp"
gdf          = gpd.read_file(fp)

gdf = gdf.assign(area=gdf.area).sort_values("area",ascending=False)






gmsh.initialize() 
gmsh.model.add("danish_waters")


factory = gmsh.model.geo
points = [] 


pnt_id  = 0
line_id = 0
curve_loop_ids = []

for poly_id, polygon in enumerate(gdf.geometry):
    
    print(f"creating p")

    print(polygon)
    x,y = polygon.exterior.xy

    # Add points
    pnt_id_start = pnt_id
    for xx,yy in zip(x[:-1],y[:-1]):
        print(f"adding point {pnt_id}")
        factory.add_point(xx, yy, 0, length_scale, pnt_id)
        pnt_id = pnt_id + 1

    # Add lines:
    line_id_start = line_id
    for i in range(pnt_id_start,pnt_id-1):
        print(f"adding point {line_id} from points {i},{i+1}")
        factory.add_line(i,i+1,line_id)
        line_id = line_id + 1
    print(f"adding point {line_id} from points {pnt_id-1},{pnt_id_start}")
    factory.add_line(pnt_id-1,pnt_id_start,line_id)

    lines = list(range(line_id_start,line_id+1))
    print(f"adding curve loop {line_id} of lines {lines}")
    factory.addCurveLoop(lines, line_id)
    curve_loop_ids.append(line_id)

    line_id = line_id + 1


    if poly_id>20:
        break

factory.addPlaneSurface(curve_loop_ids, line_id)
line_id = line_id + 1

factory.synchronize()
#gmsh.model.mesh.set_algorithm(dim=2,tag=line_id-1,val=11)
gmsh.model.mesh.generate(2)

gmsh.write("danish_waters.msh")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
