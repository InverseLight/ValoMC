function elements = findElementsInsideCircle(vmcmesh, radius, center)
% Finds elements inside a circle with a given center and radius
%
% function elements = findElementsInsideCircle(vmcmesh, radius, center)
%
% INPUT
%
%  vmcmesh:         the mesh containing the elements (described in documentation/list of structures)
%  radius:       the radius of the circle
%
% OUTPUT
%
%  center:       the center of the circle
%
   v1 = ((vmcmesh.r(vmcmesh.H(:,1),1)-center(1)).*(vmcmesh.r(vmcmesh.H(:, 1),1)- ...
                                             center(1))+ ...
         (vmcmesh.r(vmcmesh.H(:,1),2)-center(2)).*(vmcmesh.r(vmcmesh.H(:,1),2)- ...
                                             center(2)) < radius*radius);
         
   v2 = ((vmcmesh.r(vmcmesh.H(:,2),1)-center(1)).*(vmcmesh.r(vmcmesh.H(:, 2),1)- ...
                                             center(1))+ ...
         (vmcmesh.r(vmcmesh.H(:,2),2)-center(2)).*(vmcmesh.r(vmcmesh.H(:,2),2)- ...
                                             center(2)) < radius*radius);

   v3 = ((vmcmesh.r(vmcmesh.H(:,3),1)-center(1)).*(vmcmesh.r(vmcmesh.H(:, 3),1)- ...
                                             center(1))+ ...
         (vmcmesh.r(vmcmesh.H(:,3),2)-center(2)).*(vmcmesh.r(vmcmesh.H(:,3),2)- ...
                                             center(2)) < radius*radius);

   elements = find(v1==true & v2==true & v3 ==true);

end


