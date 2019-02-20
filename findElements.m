function elements = findElements(vmcmesh, querystring, varargin)
%FINDELEMENTS Finds boundary elements from the mesh
%
% USAGE:
%
%       elements = findElements(vmcmesh, querystring, varargin)
%
% DESCRIPTION:
%
%       This function can be used to find elements from the mesh.
%       A complete description is given in the homepage (see below).
%
% INPUT:
%
%       vmcmesh       - mesh structure, contains the geometry of the system
%
%       querystring, optional arguments
%
%       2D mesh (row size)
%
%          'rectangle', position(2), width (1), height (1)
%          'circle', location (2), radius (1)
%          'inverse', elements (number of elements in the selection)
%          'location', location (2)
%          'region', region_BH (2)
%
%       3D mesh (row size)
%
%          'cylinder', origin (3), direction (3), radius (1)
%          'box', origin (2), xsize (1), ysize (1), zsize (1)
%          'sphere', location (3), radius (1)
%          'inverse', elements (number of elements in the selection)
%          'location', location (3)
%          'region', region_BH (3)
%
% SEE ALSO:
%
% Detailed documentation of the function is given in
% https://inverselight.github.io/ValoMC/findingelements.html
%
% This function is provided with ValoMC

    if(size(vmcmesh.H,2) == 3)
        % 2D
        if(strcmp(querystring, 'rectangle'))
            elements = findElementsInsideRectangle(vmcmesh, varargin{2}, varargin{3}, varargin{1});
            % input is (vmcmesh, width, height, location)
        elseif(strcmp(querystring, 'circle'))
   	        elements = findElementsInsideCircle(vmcmesh, varargin{2}, varargin{1});
        elseif(strcmp(querystring, 'inverse'))
            arr = 1:size(vmcmesh.H);
            elements = setdiff(arr, [varargin{1}]);
        elseif(strcmp(querystring, 'location'))
            elements =  findElementsNearest(vmcmesh, varargin{1});
        elseif(strcmp(querystring, 'region'))
            warning('2d region not yet implemented')
        end
    elseif(size(vmcmesh.H,2) == 4)
        % 3D
        if(strcmp(querystring, 'cylinder'))
            warning('3d cylinder not yet implemented')
        elseif(strcmp(querystring, 'box'))
            warning('3d box not yet implemented')
        elseif(strcmp(querystring, 'sphere'))
            %warning('3d sphere not yet implemented')
            location = varargin{1};
            radius   = varargin{2}; 
            dist1 = vmcmesh.r(vmcmesh.H(:,1),:) - location;
            dist2 = vmcmesh.r(vmcmesh.H(:,2),:) - location;
            dist3 = vmcmesh.r(vmcmesh.H(:,3),:) - location;
            dist4 = vmcmesh.r(vmcmesh.H(:,4),:) - location;  
            norm1 = sum(dist1.^2,2);
            norm2 = sum(dist2.^2,2);
            norm3 = sum(dist3.^2,2);
            norm4 = sum(dist4.^2,2);
            r2 = radius*radius;
            elements=find(norm1  < r2 & norm2  < r2 & norm3  < r2 & norm4  < r2);
        elseif(strcmp(querystring, 'halfspace'))
            position=varargin{1};
            normal=varargin{2};
  	        elements=[];
            for i=1:size(vmcmesh.H,1)
                l1 = vmcmesh.r(vmcmesh.H(i,1),:) - position;
                l2 = vmcmesh.r(vmcmesh.H(i,2),:) - position;
                l3 = vmcmesh.r(vmcmesh.H(i,3),:) - position;
                l4 = vmcmesh.r(vmcmesh.H(i,4),:) - position;
                if(dot(l1, normal) >= 0 & dot(l2, normal) >= 0 & dot(l3, normal) >= 0 & dot(l4,normal))
                   elements = [elements; i];
                end
            end
        elseif(strcmp(querystring, 'inverse'))
            arr = 1:size(vmcmesh.H);
            elements = setdiff(arr, [varargin{1}]);
        elseif(strcmp(querystring, 'location'))
            warning('3d location not yet implemented')
        elseif(strcmp(querystring, 'region'))
            surface = varargin{1};
            indices = unique(surface(:)); % get the points that belong to surface

            minimum = min(vmcmesh.r(indices,:));
            half = (max(vmcmesh.r(indices,:)) - min(vmcmesh.r(indices,:))) / 2;
            

            DT = delaunayTriangulation((vmcmesh.r(indices,:) - minimum - half)*1.1);
            elements = [];
            for i=1:size(vmcmesh.H,1)
               if(~isnan(pointLocation(DT, vmcmesh.r(vmcmesh.H(i,1),:)-minimum-half)) && ...  
                  ~isnan(pointLocation(DT, vmcmesh.r(vmcmesh.H(i,2),:)-minimum-half)) && ...
                  ~isnan(pointLocation(DT, vmcmesh.r(vmcmesh.H(i,3),:)-minimum-half)) && ...
                  ~isnan(pointLocation(DT, vmcmesh.r(vmcmesh.H(i,4),:)-minimum-half)))
                    elements = [elements; i];
               end
            end
        end        
    else
        error('Could not recognize mesh');  
    end

end


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


function elements = findElementsInsideRectangle(vmcmesh, width, height, center)
    % Finds elements inside a rectangle
    %
    % function elements = findElementsInsideRectangle(vmcmesh, radius, center)
    %
    % INPUT
    %
    %  vmcmesh:      the vmcmesh containing the elements (described in documentation/list of structures)
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


function [segments] = findElementsNearest(vmcmesh, locations)
%
% Finds elements nearest to each position
%
% INPUT
%
%  vmcmesh:      (described in documentation/list of structures)
%  locations:    an array that contains a position vector in each row.
%
% OUTPUT
%
%  segments:     the indices of the line segments nearest to the positions
%

    avgx = (vmcmesh.r(vmcmesh.H(:,1),1) + vmcmesh.r(vmcmesh.H(:,2),1) + vmcmesh.r(vmcmesh.H(:,3),1))/3.0;
    avgy = (vmcmesh.r(vmcmesh.H(:,1),2) + vmcmesh.r(vmcmesh.H(:,2),2) + vmcmesh.r(vmcmesh.H(:,3),2))/3.0;
    pos = [avgx avgy];
    segments = zeros(size(locations,1),1);
    for ii=1:size(locations,1)
       m=(pos - locations(ii,:)) .^2;
       norms=sum(m')';
       [minvalue minindex] = min(norms);
       segments(ii) = minindex; 
    end

end

    

