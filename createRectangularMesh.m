function vmcmesh = createRectangularMesh(xsize, ysize, dh)
% Create a mesh structure with a rectangular geometry
%
% vmcmesh = createRectangularMesh(xsize, ysize, step)
%
% INPUT
%
%  xsize:       the radius of the rectangle [mm]
%  ysize:       discretisation step size [mm]
%  dh:          step size [mm]
%
% OUTPUT
%
%  vmcmesh:        structure that contains the geometry (triangles,
%                                                     boundary
%                                                     lines etc.)
%
% See also createCircularMesh
%

   vmcmesh = createGridMesh(-xsize/2+dh/2:dh:xsize/2-dh/2,-ysize/2+dh/2:dh:ysize/2-dh/2);

end

