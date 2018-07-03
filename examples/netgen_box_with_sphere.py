from netgen.csg import *
from ngsolve import Mesh


geo = CSGeometry()
# Set the mesh size on the sphere surface to 0.1
sphere = Sphere(Pnt(0,0,0),2).maxh(0.1)
# meshsize of the surface of the brick will not be finer than
# the volume mesh size
brick = OrthoBrick(Pnt(-5,-5,-5),Pnt(5,5,5))
#
#
rodbase = Cylinder(Pnt(-3, 0, 0), Pnt(0, 0, 0), 1.0);
rodplane1 = Plane(Pnt(-5,0,0), Vec(-1,0,0));
rodplane2 = Plane(Pnt(-5,0,0), Vec(1,0,0));
rod = rodbase * rodplane1 * rodplane2;

brickrod = sphere+rod;

# in the outer region we don't give a local mesh size -> global
# is used
geo.Add(brick-brickrod)
# in the volume of the sphere we set the meshsize to 0.2
geo.Add(brickrod,maxh=0.2)
# the global mesh size is set to 0.4
ngmesh = geo.GenerateMesh(maxh=0.4)

# for visualization we need a NGSolve mesh
Draw(Mesh(ngmesh))

