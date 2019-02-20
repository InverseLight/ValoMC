function elements = findBoundaries(vmcmesh, querystring, varargin)
%FINDBOUNDARIES Finds boundary elements from the mesh
%
% USAGE:
%
%       elements = findBoundaries(vmcmesh, querystring, varargin)
%
% DESCRIPTION:
%
%       This function can be used to find boundary elements from the mesh.
%       A complete description of this function is given in the homepage (see below)
%
% INPUT:
%
%       vmcmesh       - mesh structure, contains the geometry of the system
%
%       querystring, optional arguments
%
%       2D mesh (row size)
%
%          'arc', origin (2), startangle (1), endangle (1)
%          'direction', origin (2), waypoint (2), width (1)
%          'inverse', elements (number of boundary elements)
%          'location', coordinate (2) | if a third argument is given, returns nearest nodes (in the boundary) instead
%
%       3D mesh (row size)
%
%          'direction', origin (3), waypoint (3), shape ('circle','rectangle', 'arbitrary') width (1), coordinates for the shape (optional)
%          'halfspace', location (3), normal (3)
%          'inverse', elements (number of boundary elements)
%          'location', nearestlocation (3)
%
% SEE ALSO:
%
% Detailed documentation of the function is given in
%
% https://inverselight.github.io/ValoMC/findingboundaries.html
%
% This function is provided with ValoMC
    if(size(vmcmesh.H,2) == 3)
        % 2D
        if(strcmp(querystring, 'direction'))
            if(size(varargin,2) < 4)
                maxdist = ones(size(varargin{1},1),1)*1e6;
            else
                maxdist = varargin{4};
            end
            elements = findBoundariesByDirection(vmcmesh, varargin{1}, varargin{2}, varargin{3});
        elseif(strcmp(querystring, 'arc'))

            origin = varargin{1};
            startangle = varargin{2};
            endangle = varargin{3};
            % find the center of the mesh
            centerpoints1 = vmcmesh.r(vmcmesh.BH(:,1),:)-origin;
            centerpoints2 = vmcmesh.r(vmcmesh.BH(:,2),:)-origin;
            angle1=atan2(centerpoints1(:,2),centerpoints1(:,1));
            angle2=atan2(centerpoints2(:,2),centerpoints2(:,1));

            neg1=find(angle1 < 0);
            neg2=find(angle2 < 0);

            angle1(neg1) = 2*pi + angle1(neg1);
            angle2(neg2) = 2*pi + angle2(neg2);

            belements=find(angle1 >= startangle & angle1 <= endangle & angle2 >= startangle & angle2 <= endangle)

        elseif(strcmp(querystring, 'location'))
            if(size(varargin, 2) == 1)
               elements = findLineSegmentsNearest(vmcmesh, varargin{1});
            else
               elements = findNodesNearest(vmcmesh, varargin{1});
            end
        elseif(strcmp(querystring, 'rectangle'))
             
        elseif(strcmp(querystring, 'circle'))

        elseif(strcmp(querystring, 'inverse'))
            arr = 1:size(vmcmesh.BH);
            elements = setdiff(arr, [varargin{1}]);
        end
    elseif(size(vmcmesh.H,2) == 4)
        % 3D
        if(strcmp(querystring, 'direction'))
            origin=varargin{1};
            waypoint=varargin{2};
            radius=varargin{3};
            if(size(varargin, 2) < 3)
                error('Not enought arguments given') 
            end
            if(size(varargin, 2) >= 4)
                shape=varargin{4};
            else
                shape='circle';
            end
            if(size(varargin, 2) >= 5)
                points=varargin{5};
            else
            if(size(varargin,2) < 5) 
                shape = 'cylinder'; 
            else 
                shape = varargin{5};
            end
            if strcmp(shape, 'cylinder')
                dir = varargin{2} - varargin{1};
                tangent = dir/norm(dir);
                if(tangent(2) ~= tangent(1)) 
                    normal = [tangent(2) tangent(1) tangent(3)];
                else
                    normal = [tangent(1) tangent(3) tangent(2)];
                end
                normal = normal/norm(normal);
                binormal = cross(tangent, normal);
                binormal = binormal / norm(binormal);
                steps=10;
                for i=1:steps
                    points(i,:) = [radius*cos(2*pi/steps * i) radius*sin(2*pi/steps * i)];
                end  
            elseif strcmp(shape, 'rectangle') 
                    points = [-radius -radius; -radius radius; radius radius; radius -radius;]; 
            elseif strcmp(shape, 'arbitrary')
                    if(size(varargin,2) ~= 6)
                        error('Did not obtain points as input');
                    end
                    points = varargin{6};
            end
            
            transformed_points = [];
            
            for  i=1:size(points,1)
                    newpoint = origin+points(i,1)*normal + points(i,2)*binormal;
                    transformed_points = [transformed_points; newpoint];
            end
            for  i=1:size(points,1)
                    newpoint = origin+points(i,1)*normal + points(i,2)*binormal + norm(dir)*tangent;
                    transformed_points = [transformed_points; newpoint];   
            end
            
            DT = delaunayTriangulation(transformed_points);

            firstcorner = isnan(pointLocation(DT, vmcmesh.r(vmcmesh.BH(:,1),:)));
            secondcorner = isnan(pointLocation(DT, vmcmesh.r(vmcmesh.BH(:,2),:))); 
            thirdcorner = isnan(pointLocation(DT, vmcmesh.r(vmcmesh.BH(:,3),:))); 
            elements = find(firstcorner == 0 & secondcorner == 0 & thirdcorner == 0);
            end
        end
    else
            error('Could not recognize mesh');  
    end

end


function wn = windnum(line_strip, point)
   line_strip = [line_strip ; line_strip(1,:)];
   wn = 0;
   for ii=1:size(line_strip,1)-1
      if(line_strip(ii, 2) <= point(2))
         if (line_strip(ii+1, 2) > point(2))             
             if ((line_strip(ii+1, 1) - line_strip(ii, 1)) * (point(2) - line_strip(ii, 2))  - (point(1) - line_strip(ii, 1)) * (line_strip(ii+1, 2) - line_strip(ii, 2)) > 0.0)
                wn=wn+1;
             end
         end
      else
        if (line_strip(ii+1, 2) <= point(2))             
            if ((line_strip(ii+1, 1) - line_strip(ii, 1)) * (point(2) - line_strip(ii, 2))  - (point(1) - line_strip(ii, 1)) * (line_strip(ii+1, 2) - line_strip(ii, 2)) < 0.0)
                wn=wn-1;
            end
        end
      end 
   end
end



function [segments, totallength] = findBoundariesByDirection(vmcmesh, start, waypoint, width)
    segments = [];

    direction = (waypoint - start);
    prpdir = zeros(1,2);
    prpdir(1) = -direction(2);
    prpdir(2) =  direction(1);
    prpdir = prpdir / norm(prpdir);
    quad = zeros(4, 2);
    wh = width/2;
    % form a quad
    quad(1,:) = start + wh*prpdir;
    quad(2,:) = start + wh*prpdir + direction;
    quad(3,:) = start + wh*prpdir + direction - 2*wh*prpdir;
    quad(4,:) = start - wh*prpdir;
    match=1;

    for i=1:size(vmcmesh.BH, 1)
        % compute the winding number for both nodes of BH
        wn1 = windnum(quad, vmcmesh.r(vmcmesh.BH(i,1),:));
        wn2 = windnum(quad, vmcmesh.r(vmcmesh.BH(i,2),:));
        if(wn1 || wn2)
           segments = [segments; i]; 
           match = match+1;
        end
    end

    if(match == 1)
        warning('Did not find any boundaries matching the criteria.');
    end

end


function [segments] = findLineSegmentsNearest(vmcmesh, locations)
%
% Finds line segments nearest to each position
%
% INPUT
%
%  vmcmesh:         (described in documentation/list of structures)
%  locations:    an array that contains a position vector in each row.
%
% OUTPUT
%
%  segments:     the indices of the line segments nearest to the positions
%    
        avgx = (vmcmesh.r(vmcmesh.BH(:,1),1) + vmcmesh.r(vmcmesh.BH(:,2),1))/2.0;
        avgy = (vmcmesh.r(vmcmesh.BH(:,1),2) + vmcmesh.r(vmcmesh.BH(:,2),2))/2.0;
        pos = [avgx avgy];
        segments = zeros(size(locations,1),1);
        for ii=1:size(locations,1)
           m=(pos - locations(ii,:)) .^2;
           norms=sum(m')';
           [minvalue minindex] = min(norms);
           segments(ii) = minindex; 
        end
    
 end


function [nodes] = findNodesNearest(vmcmesh, locations)
        boundary_indices = unique([vmcmesh.BH(:,1) vmcmesh.BH(:,2)]);

        pos = vmcmesh.r(boundary_indices,:);
        nodes = zeros(size(locations,1),1);

        for ii=1:size(locations,1)
           m=(pos - locations(ii,:)) .^2;
           norms=sum(m')';
           [minvalue minindex] = min(norms);
           nodes(ii) = boundary_indices(minindex); 
        end
 end

function distance = distanceFromLine(p1,p2,p3)
% calculates the distance from point p3 from a line that goes trough points p1 and p2
  x0 = p3(1);
  y0 = p3(2);
  x1 = p1(1);
  y1 = p1(2);
  x2 = p2(1);
  y2 = p2(2);
  distance = abs(((y2-y1)*x0 - (x2-x1)*y0 +x2*y1- y2*x1))/norm(p2-p1);
end


