function [segments, totallength] = findBoundariesByDirection(vmcmesh, start, waypoint, width, maxdist)
% Finds line segments at a given direction from a point
%
%
% INPUT
%
%  vmcmesh:      (described in documentation/list of structures)
%  start:        starting point vector for the line definies the region
%  waypoint:     waypoint for the line that defines the region
%  width:        total width of the region
%  maxdist:      maximum distance from origin
%
% OUTPUT
%
%  segments:     the indices of the line segments within the region
%  totallength:  total lengths of the the line segments
%
%
    if(~exist('maxdist'))
        maxdist = 1e38; % big number
    end
    segments = [];
    width=width/2;
%    for j=1:size(start,1)
       match=1;
       for i=1:size(vmcmesh.BH, 1)
           % check if the two ends of the linesegment are within a
           % perpendicular distance 'width'
           if(distanceFromLine(start(:),waypoint(:),vmcmesh.r(vmcmesh.BH(i,1),:)) <= width)
               if(distanceFromLine(start(:),waypoint(:),vmcmesh.r(vmcmesh.BH(i,2),:)) <= width)
                   % if yes, check that the normal points outwards from the direction
                   direction = vmcmesh.r(vmcmesh.BH(i,2),:) - vmcmesh.r(vmcmesh.BH(i,1),:);
                   normal = [-direction(2) direction(1)];
                   
                   if((dot(normal, waypoint - start) > 0) ...
                   && (norm(vmcmesh.r(vmcmesh.BH(i,2),:) - start) < maxdist) ...
                   && (norm(vmcmesh.r(vmcmesh.BH(i,1),:) - start)  < maxdist))
                       segments = [segments; i]; 
                       match = match+1;
                   end
               end
           end
       end
%    end
    if(match == 1) 
        warning('Did not find any boundaries matching the criteria.');
    end
end


