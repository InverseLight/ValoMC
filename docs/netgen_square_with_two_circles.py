# -*- coding: utf-8 -*-

# Example for generating a simple netgen mesh

from ngsolve import *
from netgen.geom2d import SplineGeometry

geo = SplineGeometry()
p1,p12a,p12b,p2,p3,p4 = [ geo.AppendPoint(x,y) for x,y in [(-5,5),(-5,2),(-5, -2),(-5,-5),(5,-5),(5,5)] ]

geo.SetMaterial(1, "background")
geo.SetMaterial(2, "circles")

# In NetGen, regions are defined using 'left' and 'right',
# For example, leftdomain=2, rightdomain=1  means when moving from a
# from a point A to B, region 2 is always on the left side and
# region 1 is on the right.

geo.Append (["line", p1, p12a], leftdomain=1, bc="boundary1")
geo.Append (["line", p12a, p12b], leftdomain=1, bc="lightsource")
geo.Append (["line", p12b, p2], leftdomain=1, bc="boundary1")
geo.Append (["line", p2, p3], leftdomain=1, bc="boundary1")
geo.Append (["line", p3, p4], leftdomain=1, bc="boundary1")
geo.Append (["line", p4, p1], leftdomain=1, bc="boundary1")

geo.AddCircle(c=(1,1), r=1, leftdomain=2,rightdomain=1)
geo.AddCircle(c=(-2,-2), r=2, leftdomain=2,rightdomain=1)

mesh = geo.GenerateMesh(maxh=0.1)

mesh.Save("square_with_two_circles.vol")
