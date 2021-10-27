import open3d as o3d
import trimesh
import numpy as np
import pyrender
import pyglet

process = False

if process:
    pcd = o3d.io.read_point_cloud("/work/ALT/swot/swotdev/desrochesd/floodplain/run/po/output.ply")
    pcd.estimate_normals()

    # estimate radius for rolling ball
    distances = pcd.compute_nearest_neighbor_distance()
    avg_dist = np.mean(distances)
    # ~ radius = 1.5 * avg_dist   
    radius = 20 * avg_dist   

    mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(pcd, o3d.utility.DoubleVector([radius, radius * 2]))

    trimesh_obj = trimesh.Trimesh(np.asarray(mesh.vertices), np.asarray(mesh.triangles), vertex_normals=np.asarray(mesh.vertex_normals))
    trimesh_obj.export('/work/ALT/swot/swotdev/desrochesd/floodplain/run/po/stuff.stl')





fuze_trimesh = trimesh.load('/work/ALT/swot/swotdev/desrochesd/floodplain/run/po/test_pyvista30_smooth.stl')
# ~ fuze_trimesh = trimesh.load('/work/ALT/swot/swotdev/desrochesd/trimesh/models/fuze.obj')

                              
                                      
mesh = pyrender.Mesh.from_trimesh(fuze_trimesh, smooth=False)
scene = pyrender.Scene()
scene.add(mesh)
pyrender.Viewer(scene)
