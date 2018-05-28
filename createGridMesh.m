function output = createGridMesh(xvec, yvec, zvec)
%CREATEGRIDMESH Creates a triangular/tetrahedral finite element mesh that consists of pixels/voxels
%
%
% DESCRIPTION:
%
%       This function can be used to create a mesh compatible
%       for pixel or voxel based input values. The structure of the mesh
%       is given in the figure below and has the same idea as the
%       native MATLAB function 'meshgrid'
%
%
%       o-------o--------o
%       |       |        |
%       |   x   |   x    |
%       |       |        |
%       o-------o--------o
%       |       |        |
%       |   x   |   x    |
%       |       |        |
%       o-------o--------o
%
%       o = mesh cooridinates
%       x = grid point (locattion given by xvec, yvec)
%
%       The mesh can be either 2d  (xvec, yvec given) or 3d.
%       In 2d, the triangular mesh is constructed in the following fashion   
%
%       o--b3---o--b4---o
%       | .  t6 | .  t8 |
%       b2  .   |   .   b5
%       | t2  . |t4   . |
%       o-------o-------o
%       | .  t5 | .  t7 |
%       b1  .   |   .   b6
%       | t1  . |t3   . |
%       o---b8--o---b7--o
%
%       t_i = triangle
%       b_i = boundary segment
%
%       The idea for the 3d mesh is the same. xvec, yvec and zvec
%       depict the center location of each cube. Each cube contains 6
%       tetrahedrons. The first tetrahedon is in the first cube, the
%       second tetrahedron is in the second cube, in ascending
%       y-x-z. The second tetrahedron of the first cube is the
%       size(xvec)*size(yvec)*size(xvec)+1:th tetrahedron, the third
%       tetrahedron of the first cube cube is the
%       2*size(xvec)*size(yvec)*size(xvec)+1:th tetrahedron and so on
%
% USAGE:
%       vmcmesh = createGridMesh(xvec, yvec);
%       vmcmesh = createGridMesh(xvec, yvec);
%
% INPUTS:
%       xvec       - x coordinates for the center location of each pixel/cube
%       yvec       - y coordinates for the center location of each pixel/cube
%
% OPTIONAL INPUTS:
%       zvec       - z coordinates for the center location of each pixel/cube 


    if(nargin==2)
   
        output.r = [];
        output.H = [];
        output.BH = [];
	
     	dx = abs(xvec(2)-xvec(1));
    	dy = abs(yvec(2)-yvec(1));

        gridvecx = (xvec - dx/2);
        gridvecy = (yvec - dy/2);
        
        if(iscolumn(gridvecx))
            gridvecx = [gridvecx; gridvecx(end)+dx];
            gridvecy = [gridvecy; gridvecy(end)+dy];
        else 
            gridvecx = [gridvecx gridvecx(end)+dx];
            gridvecy = [gridvecy gridvecy(end)+dy];            
        end
        
        [x,y] = meshgrid(gridvecx, gridvecy);
        output.r = [x(:) y(:)];            
        
        num_x_voxels = length(xvec);
        num_y_voxels = length(yvec);
    
        k=1;
        % k------k+1
        % | t1  . | 
        % |   .   | 
        % | .  t2 |
        %k+xwidth-k+xwidth+1    

        % Build triangles
        %
        %  3-------6--------9
        %  |       |        |
        %  |   x   |   x    |
        %  |       |        |
        %  2-------5--------8
        %  |       |        |
        %  |   x   |   x    |
        %  |       |        |
        %  1-------4--------7    
        %
        %  num_y_voxels = 2;
        %  ysize  = 3;
        

        % TODO: for loops can be optimized away
        ysize=num_y_voxels+1;
        xsize=num_x_voxels+1;
        
        output.H = zeros(num_x_voxels*num_y_voxels*2,3);
        n=1;
        for i = 1:num_x_voxels
            for j = 1:num_y_voxels
                output.H(n,:) = [k,k+1,k+ysize];
                k=k+1;
                n=n+1;
            end
            k=k+1;
        end
        k=1;
        for i = 1:num_x_voxels
            for j = 1:num_y_voxels
                output.H(n,:) = [k+1,k+ysize+1,k+ysize];
                n=n+1;
                k=k+1;
            end
            k=k+1;
        end

        output.BH = [];

        % Build boundary
        %                           
        %  TL               TR       
        %  4-------8--------12     
        %  |       |        |          
        %  |   x   |   x    |            
        %  |       |        |       
        %  3-------7--------11
        %  |       |        |
        %  |   x   |   x    |
        %  |       |        |
        %  2-------6--------10
        %  |       |        |
        %  |   x   |   x    |
        %  |       |        |
        %  1-------5--------9       
        %  BL               BR
        %
        % xsize = 3
        % ysize = 4
        
        % BL to TL
        for j = 1:num_y_voxels
            output.BH = [output.BH; [j,j+1]];
        end

        % TL to BR
        for j = 1:num_x_voxels
            output.BH = [output.BH; [ysize*j,ysize*(j+1)]];
        end

        % TR to BR
        for j = 1:num_y_voxels
            output.BH = [output.BH; [ysize*xsize - (j-1), ysize*xsize - (j) ]];
        end

        % BR to TL
        for j = 1:num_x_voxels
            output.BH = [output.BH; [ysize*(xsize-1)+1-ysize*(j-1),ysize*(xsize-1)+1-ysize*(j)]];
        end
    else
        output.r = [];
        output.H = [];
        output.BH = [];
	
     	dx = abs(xvec(2)-xvec(1));
    	dy = abs(yvec(2)-yvec(1));
    	dz = abs(zvec(2)-zvec(1));

        gridvecx = (xvec - dx/2);
        gridvecy = (yvec - dy/2);
        gridvecz = (zvec - dz/2);
        
        if(iscolumn(gridvecx))
            gridvecx = [gridvecx; gridvecx(end)+dx];
            gridvecy = [gridvecy; gridvecy(end)+dy];
            gridvecz = [gridvecz; gridvecz(end)+dz];
        else 
            gridvecx = [gridvecx gridvecx(end)+dx];
            gridvecy = [gridvecy gridvecy(end)+dy];            
            gridvecz = [gridvecz gridvecz(end)+dz];
        end
        
        [x,y,z] = meshgrid(gridvecx, gridvecy,gridvecz);
        output.r = [x(:) y(:) z(:)];     
        
        points = [0 0 0; 0 1 0; 1 0 0; 1 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 1];
        connectivity =  [5 1 2 3; 6 5 2 3; 6 7 5 3; 6 4 7 3; 6 2 4 3; 6 8 7 4];        
        
        ijk_to_index = @(i,j,k) ((k-1)*(size(x,2)*size(x,1))+(i-1)*size(x,1)+j);
        num_x_voxels=size(x,2)-1; 
        num_y_voxels=size(x,1)-1; 
        num_z_voxels=size(x,3)-1; 
        output.H=zeros(num_x_voxels*num_y_voxels*num_z_voxels*6,4);

        % Build the topology matrix
        counter=1;
        for t=1:6
            for k=1:num_z_voxels
                for i=1:num_x_voxels
                    for j=1:num_y_voxels
                       offset = connectivity(t,:);
                       ijkoffset = points(offset,:);
                       output.H(counter,:) = [ijk_to_index(i+ijkoffset(1,1),j+ijkoffset(1,2),k+ijkoffset(1,3))
                                              ijk_to_index(i+ijkoffset(2,1),j+ijkoffset(2,2),k+ijkoffset(2,3))
                                              ijk_to_index(i+ijkoffset(3,1),j+ijkoffset(3,2),k+ijkoffset(3,3))
                                              ijk_to_index(i+ijkoffset(4,1),j+ijkoffset(4,2),k+ijkoffset(4,3))];
                       counter=counter+1;
                    end   
                end
            end
        end
        
        % Build the boundary topology
        counter = 1;
        output.BH=zeros(num_x_voxels*num_y_voxels*2*2+num_y_voxels*num_z_voxels*2*2+num_z_voxels*num_x_voxels*2*2,3);

        % XY plane
        for i=1:num_x_voxels
            for j=1:num_y_voxels
               k=1;
               output.BH(counter, :) = [ijk_to_index(i,j,k) ijk_to_index(i+1,j,k) ijk_to_index(i,j+1,k)];
               counter = counter + 1;
               output.BH(counter, :) = [ijk_to_index(i+1,j+1,k) ijk_to_index(i,j+1,k) ijk_to_index(i+1,j,k)];
               counter = counter + 1;
               k=num_z_voxels+1;      
               output.BH(counter, :) = [ijk_to_index(i,j,k) ijk_to_index(i+1,j,k) ijk_to_index(i,j+1,k)];
               counter = counter + 1;
               output.BH(counter, :) = [ijk_to_index(i+1,j+1,k) ijk_to_index(i,j+1,k) ijk_to_index(i+1,j,k)];
               counter = counter + 1;
            end
        end

        % YZ plane
        for j=1:num_y_voxels
            for k=1:num_z_voxels
               i=1;
               output.BH(counter, :) = [ijk_to_index(i,j,k) ijk_to_index(i,j+1,k) ijk_to_index(i,j,k+1)];
               counter = counter + 1;
               output.BH(counter, :) = [ijk_to_index(i,j+1,k+1) ijk_to_index(i,j,k+1) ijk_to_index(i,j+1,k)];
               counter = counter + 1;
               i=num_x_voxels+1;
               output.BH(counter, :) = [ijk_to_index(i,j,k) ijk_to_index(i,j+1,k) ijk_to_index(i,j,k+1)];
               counter = counter + 1;
               output.BH(counter, :) = [ijk_to_index(i,j+1,k+1) ijk_to_index(i,j,k+1) ijk_to_index(i,j+1,k)];
               counter = counter + 1;
            end
        end

        % ZX plane
        for k=1:num_z_voxels
            for i=1:num_x_voxels
               j=1;
               output.BH(counter, :) = [ijk_to_index(i,j,k) ijk_to_index(i,j,k+1) ijk_to_index(i+1,j,k)];
               counter = counter + 1;
               output.BH(counter, :) = [ijk_to_index(i+1,j,k+1) ijk_to_index(i+1,j,k) ijk_to_index(i,j,k+1)];
               counter = counter + 1;
               j=num_y_voxels+1;
               output.BH(counter, :) = [ijk_to_index(i,j,k) ijk_to_index(i,j,k+1) ijk_to_index(i+1,j,k)];
               counter = counter + 1;
               output.BH(counter, :) = [ijk_to_index(i+1,j,k+1) ijk_to_index(i+1,j,k) ijk_to_index(i,j,k+1)];
               counter = counter + 1;
            end
        end


end
