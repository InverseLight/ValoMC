from netgen.csg import *
from ngsolve import Mesh


geo = CSGeometry()
# Set the mesh size on the sphere surface to 0.1
sphere = Sphere(Pnt(0,0,0),2).maxh(0.4)

brick = OrthoBrick(Pnt(-5,-5,-5),Pnt(5,5,5))


rodbase = Cylinder(Pnt(-3, 0, 0), Pnt(0, 0, 0), 1.0);
rodplane1 = Plane(Pnt(-5,0,0), Vec(-1,0,0));
rodplane2 = Plane(Pnt(-5,0,0), Vec(1,0,0));
rod = rodbase * rodplane1 * rodplane2;

brickrod = sphere+rod;

geo.Add(brick-brickrod)

geo.Add(brickrod,maxh=0.4)

ngmesh = geo.GenerateMesh(maxh=0.8)

ngmesh.Save("sphere_in_box.vol")

