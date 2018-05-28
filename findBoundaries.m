function elements = findBoundaries(vmcmesh, querystring, varargin)
% Finds elements
%
% See "Finding elements" in documentation
    if(size(vmcmesh.H,2) == 3)
        % 2D
        if(strcmp(querystring, 'direction'))
            if(size(varargin,2) < 4)
                maxdist = ones(size(varargin{1},1),1)*1e6;
            else
                maxdist = varargin{4};
            end
            elements = findBoundariesByDirection(vmcmesh, varargin{1}, varargin{2}, varargin{3}, maxdist);
        elseif(strcmp(querystring, 'location'))
            elements = findLineSegmentsNearest(vmcmesh, varargin{1});
        elseif(strcmp(querystring, 'rectangle'))

        elseif(strcmp(querystring, 'circle'))

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
    else
            error('Could not recognize mesh');  
    end

end


