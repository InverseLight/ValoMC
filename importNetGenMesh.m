function [vmcmesh regions region_names boundaries boundary_names] = importNetGenMesh(filename_vol)
% Import a NetGen mesh
%
% function [mesh regions region_names boundaries boundary_names] = importNetGenMesh(filename_vol)
%
% INPUT
%
%  filename:            filename of a file that is in Netgen's native 'vol' -format
%
% OUTPUT
%
%  mesh:		mesh structure (see the documentation for the structures in the toolbox doc)
%  regions: 		the regions in the vol file as cell arrays. Each element contains a vector
%                       that holds the indices of a region in the 'medium' structure.
%  region_names:        names of the regions as a cell array. Indexing is the same as in 'regions',
%                       i.e. region_names(i) is the name of the region indicated by regions(i)
%  boundaries: 		same as regions, but for boundaries
%  boundary_names:	same as region_names, but for boundaries
%
%
%  EXAMPLE:
%
% [mesh regions region_names boundaries boundary_names] = importNetGenMesh('mymesh.vol');
% indices_for_background = cell2mat(regions(1));
% indices_for_circles = cell2mat(regions(2));
% indices_for_lightsource = cell2mat(boundaries(2));
% vmcmedium.refractive_index(indices_for_circles)=1.6;
% vmcboundary.lightsource(indices_for_lightsource) = {'cosinic'};

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

    else if(dimensionality == 3)
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
    else
        error('Could not read dimensionality!');
    end
end 
