function [vmcmesh regions region_names boundaries boundary_names] = importNetGenMesh(filename_vol, preserve_boundary)
% IMPORTNETGENMESH Imports NetGen .vol files
%
% USAGE
%
%  function [mesh regions region_names boundaries boundary_names] = importNetGenMesh(filename_vol)
%
% INPUT
%
%  filename:            filename ('.vol', file must be in ASCII format)
%
% OPTIONAL INPUT
%
%  preserve_boundary:   Set to false to build a new boundary (BH) for the mesh. 
%                       This ignores the boundary definition in the file. 
%                       The new boundary is build so that there no boundary elements between elements. 
%
% OUTPUT
%
%  mesh:                a structure that contains coordinates (r), elements (H, indices to r) 
%                       and boundary elements (BH, indices to r) of the mesh
%
%  regions: 		    the regions (indices to H) in the vol file as cell array
%
%  region_names:        names of the regions as a cell array. Indexing is the same as in 'regions'.
%
%  boundaries: 		    same as regions, but for boundaries (indices to BH)
%
%  boundary_names:	    same as region_names, but for boundaries
%
%
%  EXAMPLE:
%
%   [mesh regions region_names boundaries boundary_names] = importNetGenMesh('mymesh.vol');
%   indices_for_background = cell2mat(regions(1));
%   indices_for_circles = cell2mat(regions(2));
%   indices_for_lightsource = cell2mat(boundaries(2));
%   vmcmedium.refractive_index(indices_for_circles)=1.6;
%   vmcboundary.lightsource(indices_for_lightsource) = {'cosinic'};
%
%   This function is provided with ValoMC
    

    fid = -1;
    errmsg = '';
    [fid,errmsg] = fopen(filename_vol);
    if(fid < 0)
        disp(errmsg);
        return
    end
    surface_elements = [];
    volume_elements = [];
    materials = [];
    bcnames = [];
    points = [];
    edgesegments = [];
    points = [];
    dimensionality = [];
    tline = ' ';
    while ischar(tline)
        % read the file line by line and scan for arrays
        tline = fgetl(fid);
        if(isempty(dimensionality))
            dimensionality = readNetgenValue('dimension', fid,tline);
        end
        if(isempty(surface_elements))
            surface_elements = readNetgenEntry('surfaceelements', fid,tline);
        end 
        if(isempty(surface_elements))
            surface_elements = readNetgenEntry('surfaceelementsgi', fid,tline);
        end
        if(isempty(volume_elements))
            volume_elements = readNetgenEntry('volumeelements', fid,tline);
        end 
        if(isempty(volume_elements))
            volume_elements = readNetgenEntry('volumeelementsgi', fid,tline);
        end 
        if(isempty(materials))
            materials = readNetgenEntry('materials', fid,tline);
        end
        if(isempty(bcnames))
            bcnames = readNetgenEntry('bcnames', fid,tline);
        end
        if(isempty(edgesegments))
            edgesegments = readNetgenEntry('edgesegmentsgi2', fid,tline);
        end
        if(isempty(edgesegments))
            edgesegments = readNetgenEntry('edgesegments', fid,tline);
        end
        if(isempty(points))
            points = readNetgenEntry('points', fid,tline);    
        end
    end  
    fclose(fid);    
    if(dimensionality == 2)  
        vmcmesh.H = surface_elements(:,6:8);
        vmcmesh.r = points(:,1:2);
        vmcmesh.BH =edgesegments(:,3:4);
        region_names = materials;
        boundary_names = bcnames;
        num_boundaries = max(edgesegments(:,1));
        num_domains = max(surface_elements(:,2));

        for i = 1:num_boundaries
            boundaries{i} = find(edgesegments(:,1)==i);
        end 
        for i = 1:num_domains
            regions{i} = find(surface_elements(:,2)==i);
        end 

        if(exist('preserve_boundary'))
            if(~preserve_boundary)              
               vmcmesh.BH = sort(createBH(vmcmesh.H,createHN(vmcmesh.H)),2); 
               for i = 1:num_boundaries 
                  cur_elements=sort(cell2mat(boundaries{i}),2);
                  [~,~,new_indices] = intersect(cur_elements,vmcmesh.BH,'rows');
                  boundaries{i} = new_indices;
               end
            end
        end    
    
    elseif(dimensionality == 3)
        vmcmesh.H = volume_elements(:,3:6);
        vmcmesh.r = points(:,1:3);
        vmcmesh.BH = surface_elements(:,6:8);
        region_names = materials;
        boundary_names = bcnames;
        num_boundaries = max(surface_elements(:,2));
        num_domains = max(volume_elements(:,1));        
        for i = 1:num_boundaries
            boundaries{i} = find(surface_elements(:,2)==i);
        end         
        for i = 1:num_domains    
            regions{i} = find(volume_elements(:,1)==i);
        end         
        if(exist('preserve_boundary'))
            if(~preserve_boundary)              
               vmcmesh.BH = sort(createBH(vmcmesh.H),2); 
               for i = 1:num_boundaries 
                  cur_elements=sort(cell2mat(boundaries{i}),2);
                  [~,~,new_indices] = intersect(cur_elements,vmcmesh.BH,'rows');
                  boundaries{i} = new_indices;
               end
            end
        end    
    else
        error('Could not read dimensionality!');
    end



end


function [arr_out] = readNetgenEntry(entryname, fid, tline)
    arr_out = [];
    value = [];
    if(strcmp(tline, entryname))
        number_of_entries = str2num(fgetl(fid));
        for i=1:number_of_entries
            curline = fgetl(fid);
            linevector = str2num(curline);
            if(isempty(linevector))
                [value] = strsplit(curline);
                arr_out{uint32(str2num(value{1}))} = value{2};
            else
                arr_out = [arr_out; linevector];
            end
        end
        if(not(ischar(tline)))
            error('Could not read file');
        end
    end
end


function [value] = readNetgenValue(entryname, fid, tline)
% Auxiliary function for importNetGenmesh
%
%  This function is used to extract scalars from the vol files
%  (like the dimensionality).
%
    value = [];
    if(strcmp(tline, entryname))
        value = str2num(fgetl(fid));
    end
end


