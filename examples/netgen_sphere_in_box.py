from netgen.csg import *
from ngsolve import Mesh


geo = CSGeometry()

# Create a box as background
brick = OrthoBrick(Pnt(-5,-5,-5),Pnt(5,5,5))
brick.mat("background");
# Add a sphere into the middle
sphere = Sphere(Pnt(0,0,0),2).maxh(0.4)
sphere.mat("sphere")

# Remesh a face of the cube so that it contains
# a cylindrical inclusion

infcylinder = Cylinder(Pnt(-5, 0, 0), Pnt(5, 0, 0), 1.0); # infinite cylinder along x = -5 ... 5
# Create two planes to make that cut the cylinder to make it finite
plane1 = infcylinder*Plane (Pnt(-5,0,0), Vec(-1,0,0) );
plane2 = infcylinder*Plane (Pnt(-4,0,0), Vec(1,0,0) );
cylinder = (plane1*plane2*infcylinder).bc("cylinder");

geo.Add(cylinder+brick-sphere);
geo.Add(sphere);

ngmesh = geo.GenerateMesh(maxh=0.4)
ngmesh.Save("sphere_in_box.vol")

