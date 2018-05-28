function exportX3D(filename, vertices, mesh, value, colormap)
%EXPORTX3D Export the mesh into X3D format.
%
% DESCRIPTION:
%
%       Export the mesh to X3D format for use in e.g. 3d modeling
%       software.  The mesh has to be a three dimensional, i.e. to
%       consists of tetrahedrons. Each tetrahedron will be stored as 4
%       triangles.  The triangles will be colored according to
%       'colormap' and 'value'. The length of 'value' is equal to the
%       number of tetrahedrons
%
% USAGE:
%
%       exportX3D('mouse.x3d', vmcmesh.r,vmcmesh.H, solution.element_fluence, @parula)
%       exportX3D('mouse.xhtml', vmcmesh.r, vmcmesh.H, solution.element_fluence, @jet)
%
% INPUTS:
%
%       filename       - see usage example
%       vertices       - vector for node coordinates
%       mesh           - tetrahederal mesh
%       value          - value for each tetrahedron
%       colormap       - colormap function for determining the
%                        color of each triangle from 'value'
%
%

% determine if xhtml header is needed

[~,~,ext] = fileparts(filename);

is_xhtml = false;

if(strcmp(ext,'.xhtml')) 
    is_xhtml=true;
end
% create an index array for triangles

n_tetrahedrons=size(mesh,1); 
n_triangles = size(mesh,1)*4;



triangles(1:4:n_triangles,1) = mesh(1:n_triangles/4,1); % 1st triangle 
triangles(1:4:n_triangles,2) = mesh(1:n_triangles/4,2);
triangles(1:4:n_triangles,3) = mesh(1:n_triangles/4,3);

triangles(2:4:n_triangles,1) = mesh(1:n_triangles/4,1); % 2st triangle 
triangles(2:4:n_triangles,2) = mesh(1:n_triangles/4,2);
triangles(2:4:n_triangles,3) = mesh(1:n_triangles/4,4);

triangles(3:4:n_triangles,1) = mesh(1:n_triangles/4,1); % 3st triangle 
triangles(3:4:n_triangles,2) = mesh(1:n_triangles/4,3);
triangles(3:4:n_triangles,3) = mesh(1:n_triangles/4,4);

triangles(4:4:n_triangles,1) = mesh(1:n_triangles/4,2); % 4st triangle 
triangles(4:4:n_triangles,2) = mesh(1:n_triangles/4,3);
triangles(4:4:n_triangles,3) = mesh(1:n_triangles/4,4);





triangle_center = (vertices(triangles(1:n_triangles,1),:)+ ...
                            vertices(triangles(1:n_triangles,2),:)+ ...
                            vertices(triangles(1:n_triangles,3),:))/3;

% create colors for triangles

% open file

folder=which('exportX3D');

if(size(folder,1)==0)
    warning('Could not find exportX3D from the search path.');
    return 
end

[filepath,~,~] = fileparts(folder);

if(is_xhtml) 
    textbase = fileread([filepath '/basis.xhtml']);    
else
    textbase = fileread([filepath '/basis.x3d']);        
end

basisstring = string(textbase);

% create a colormap

color_matrix = colormap(256);

% convert coordinates and indices to characters

r_vector(1:3:size(vertices,1)*3) = vertices(:,1);
r_vector(2:3:size(vertices,1)*3) = vertices(:,2);
r_vector(3:3:size(vertices,1)*3) = vertices(:,3);





tetrahedron_center = (vertices(mesh(:,1),:)+...
                                  vertices(mesh(:,2),:)+...
                                  vertices(mesh(:,3),:)+...
                                  vertices(mesh(:,4),:))/4;

                             
                             
for ii = 1:n_tetrahedrons
% make sure the normal of each triangle points outward from the 
% center of the tetrahedron
   triangle_index = (ii-1)*4 + 1;   
   
   tetrahedron_center2 = (vertices(mesh(ii,1),:)+...
                                       vertices(mesh(ii,2),:)+...
                                       vertices(mesh(ii,3),:)+...
                                       vertices(mesh(ii,4),:))/4;

 %  v2=vertices(triangles(triangle_index,2),:) - vertices(triangles(triangle_index,1),:); 
  % v3=vertices(triangles(triangle_index,3),:) - vertices(triangles(triangle_index,1),:); 
  
   triangle_normal1=cross(vertices(triangles(triangle_index,2),:) - vertices(triangles(triangle_index,1),:), vertices(triangles(triangle_index,3),:) - vertices(triangles(triangle_index,1),:));
   triangle_normal2=cross(vertices(triangles(triangle_index+1,2),:) - vertices(triangles(triangle_index+1,1),:), vertices(triangles(triangle_index+1,3),:) - vertices(triangles(triangle_index+1,1),:));
   triangle_normal3=cross(vertices(triangles(triangle_index+2,2),:) - vertices(triangles(triangle_index+2,1),:), vertices(triangles(triangle_index+2,3),:) - vertices(triangles(triangle_index+2,1),:));
   triangle_normal4=cross(vertices(triangles(triangle_index+3,2),:) - vertices(triangles(triangle_index+3,1),:), vertices(triangles(triangle_index+3,3),:) - vertices(triangles(triangle_index+3,1),:));

                                   
   t1=tetrahedron_center(ii,:) - triangle_center(triangle_index,:);
   t2=tetrahedron_center(ii,:) - triangle_center(triangle_index+1,:);
   t3=tetrahedron_center(ii,:) - triangle_center(triangle_index+2,:);
   t4=tetrahedron_center(ii,:) - triangle_center(triangle_index+3,:);

   
   if dot(triangle_normal1, t1) > 0
       swp=triangles(triangle_index,3);
       triangles(triangle_index,3) = triangles(triangle_index,2);
       triangles(triangle_index,2)=swp;
   end

   if dot(triangle_normal2, t2) > 0
       swp=triangles(triangle_index+1,3);
       triangles(triangle_index+1,3) = triangles(triangle_index+1,2);
       triangles(triangle_index+1,2)=swp;
   end

   if dot(triangle_normal3, t3) > 0
       swp=triangles(triangle_index+2,3);
       triangles(triangle_index+2,3) = triangles(triangle_index+2,2);
       triangles(triangle_index+2,2)=swp;
   end
   
   if dot(triangle_normal4, t4) > 0
       swp=triangles(triangle_index+3,3);
       triangles(triangle_index+3,3) = triangles(triangle_index+3,2);
       triangles(triangle_index+3,2)=swp;
   end
   
   
   triangle_normal1=cross(vertices(triangles(triangle_index,2),:) - vertices(triangles(triangle_index,1),:), vertices(triangles(triangle_index,3),:) - vertices(triangles(triangle_index,1),:));
   triangle_normal2=cross(vertices(triangles(triangle_index+1,2),:) - vertices(triangles(triangle_index+1,1),:), vertices(triangles(triangle_index+1,3),:) - vertices(triangles(triangle_index+1,1),:));
   triangle_normal3=cross(vertices(triangles(triangle_index+2,2),:) - vertices(triangles(triangle_index+2,1),:), vertices(triangles(triangle_index+2,3),:) - vertices(triangles(triangle_index+2,1),:));
   triangle_normal4=cross(vertices(triangles(triangle_index+3,2),:) - vertices(triangles(triangle_index+3,1),:), vertices(triangles(triangle_index+3,3),:) - vertices(triangles(triangle_index+3,1),:));

                                   
   t1=tetrahedron_center(ii,:) - triangle_center(triangle_index,:);
   t2=tetrahedron_center(ii,:) - triangle_center(triangle_index+1,:);
   t3=tetrahedron_center(ii,:) - triangle_center(triangle_index+2,:);
   t4=tetrahedron_center(ii,:) - triangle_center(triangle_index+3,:);
   
end


H_vector(1:4:size(triangles,1)*4) = int32(triangles(:,1));
H_vector(2:4:size(triangles,1)*4) = int32(triangles(:,2));
H_vector(3:4:size(triangles,1)*4) = int32(triangles(:,3));
H_vector(4:4:size(triangles,1)*4) = int32(0);

color_vector(1:3:size(color_matrix,1)*3) = color_matrix(:,1);
color_vector(2:3:size(color_matrix,1)*3) = color_matrix(:,2);
color_vector(3:3:size(color_matrix,1)*3) = color_matrix(:,3);

normalised_value = int64(255*(value - min(value))/(max(value)-min(value)));

normalised_value2(1:4:n_triangles) =  normalised_value;
normalised_value2(2:4:n_triangles) =  normalised_value;
normalised_value2(3:4:n_triangles) =  normalised_value;
normalised_value2(4:4:n_triangles) =  normalised_value;

vertice_characters = sprintf('%6f %6f %6f, ', r_vector);
H_characters = sprintf('%i ', H_vector-1);
value_characters = sprintf('%i ', normalised_value2);
color_characters = sprintf('%6f %6f %6f, ', color_vector);

basisstring = strrep(basisstring,'$VERTICES',vertice_characters);
basisstring = strrep(basisstring,'$COORDINDEX', H_characters);
basisstring = strrep(basisstring,'$COLORINDEX', value_characters);
basisstring = strrep(basisstring,'$COLORS', color_characters);

% print the file

fileID = fopen(filename,'w');
fprintf(fileID,'%s',basisstring);
fclose(fileID);

end


