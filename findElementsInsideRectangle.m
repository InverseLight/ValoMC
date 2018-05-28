function elements = findElementsInsideRectangle(vmcmesh, width, height, center)
% Finds elements inside a rectangle
%
% function elements = findElementsInsideRectangle(vmcmesh, radius, center)
%
% INPUT
%
%  vmcmesh:         the vmcmesh containing the elements (described in documentation/list of structures)
%  width:        the width of the circle
%  height:       the height of the circle
%  center:       the center of the circle
%
% OUTPUT
%
%  elements:     the elements of the rectangle

   maxx = center(1)+width/2;
   minx = center(1)-width/2;
   maxy = center(2)+height/2;
   miny = center(2)-height/2;
 
   elements = find(vmcmesh.r(vmcmesh.H(:,1),1) <= maxx & vmcmesh.r(vmcmesh.H(:,1),1) >= minx & ...
                   vmcmesh.r(vmcmesh.H(:,2),1) <= maxx & vmcmesh.r(vmcmesh.H(:,2),1) >= minx & ...
                   vmcmesh.r(vmcmesh.H(:,3),1) <= maxx & vmcmesh.r(vmcmesh.H(:,3),1) >= minx & ...
     	           vmcmesh.r(vmcmesh.H(:,1),2) <= maxy & vmcmesh.r(vmcmesh.H(:,1),2) >= miny & ...
                   vmcmesh.r(vmcmesh.H(:,2),2) <= maxy & vmcmesh.r(vmcmesh.H(:,2),2) >= miny & ...
                   vmcmesh.r(vmcmesh.H(:,3),2) <= maxy & vmcmesh.r(vmcmesh.H(:,3),2) >= miny);

end


