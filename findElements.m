function elements = findElements(vmcmesh, querystring, varargin)
% Finds elements
%
% See https://inverselight.github.io/ValoMC/findingelements.html
    if(size(vmcmesh.H,2) == 3)
        % 2D
        if(strcmp(querystring, 'rectangle'))
            elements = findElementsInsideRectangle(vmcmesh, varargin{2}, varargin{3}, varargin{1});
            % input is (vmcmesh, width, height, location)
        elseif(strcmp(querystring, 'circle'))
   	        elements = findElementsInsideCircle(vmcmesh, varargin{1}, varargin{2});
        elseif(strcmp(querystring, 'inverse'))
            arr = 1:size(vmcmesh.H);
            elements = setdiff(arr, [varargin{1}]);
        elseif(strcmp(querystring, 'location'))
            warning('2d location not yet implemented')
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


