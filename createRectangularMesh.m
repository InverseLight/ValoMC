function vmcmesh = createRectangularMesh(xsize, ysize, dh)
% Creates a mesh structure with a rectangular geometry
%
% DESCRIPTION:
%
%       This can function can be used to create a rectangular mesh (2D)
%
% USAGE:
%
%       vmcmesh = createRectangularMesh(xsize, ysize, step)
%       
% INPUT:
%
%       xsize        - rectangle width [mm]
%       ysize        - rectangle height [mm]
%       dh           - width and height of each triangle in the mesh
%
% OUTPUT:
%
%       vmcmesh      - mesh structure, contains the geometry of the system
%
% SEE ALSO:
% 
%                      createCircularMesh, createGridMesh
%
% This function is provided with ValoMC

   vmcmesh = createGridMesh(-xsize/2+dh/2:dh:xsize/2-dh/2,-ysize/2+dh/2:dh:ysize/2-dh/2);

end

