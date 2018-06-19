function elements = findElementsInsideCircle(vmcmesh, radius, center)
%FINDELEMENTSINSIDECIRCLE Returns indices to elements within a given radius from a location
%
% DESCRIPTION:
%       Returns elements within given radius from a location
%
% USAGE:
%       elements = findElementsInsideCircle(vmcmesh, center, radius)
%
% INPUTS:
%       vmcmesh     - https://inverselight.github.io/ValoMC/structures.html
%       radius      - radius of the circle
%       center      - location vector of the circle
%                     
% OUTPUTS:
%       elements    - elements within the circle

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


