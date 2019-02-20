function [vmcboundary, vmcmesh_out] = createBoundary(vmcmesh, vmcmedium_in, vmcboundary_in)
% Creates a boundary structure for a mesh
%
%  USAGE:
%
%       vmcboundary = createBoundary(vmcmesh)
%
%           Create a boundary structure with a minimal number of fields
%           in order to run a simulation, i.e. a lightsource field whose length is
%           equal to the number of boundary elements
%
%        vmcboundary = createBoundary(vmcmesh, vmcmedium)
%
%           In addition, create exterior_refactive_index so that
%           photons experience no change in refractive index as they
%           leave the medium
%
%        vmcboundary_out = createBoundary(vmcmesh, vmcmedium, vmcboundary)
%
%           Repeat the entries in vmcboundary so that the lengths
%           of the arrays are equal to the number of boundary elements in vmcmesh.
%
%        [vmcboundary, vmcmesh_out] = createBoundary(vmcmesh, ...)
%
%           In addition, create a boundary element topology (BH) in vmcmesh_out
%
% INPUT:
%
%       vmcmesh      - mesh structure, contains the geometry of the system
%
% OPTIONAL INPUT:
%
%       vmcmedium    - medium structure, used to set e. g. optical coefficients
%       vmcboundary  - boundary structure, used to set e. g. light sources
%
% OUTPUT:
%
%       vmcboundary  - boundary structure
%
% OPTIONAL OUTPUT:
%
%       vmcmesh_out  - vmcmesh complemented with boundary element structure
%
% SEE ALSO:
%
%       https://inverselight.github.io/ValoMC/structures.html
%
% This function is provided with ValoMC

    if(nargout == 2)
        vmcmesh_out = vmcmesh;
        if(~isfield(vmcmesh, 'BH'))
           vmcmesh_out.BH = createBH(vmcmesh.H, createHN(vmcmesh.H));
        end
    end

    if(exist('vmcboundary_in'))
       if(~isfield(vmcboundary_in, 'lightsource'))
           vmcboundary.lightsource = cell(size(vmcmesh.BH,1), 1);
       else
           vmcboundary.lightsource = vmcboundary_in.lightsource;
       end
    else
       vmcboundary = struct();
       vmcboundary.lightsource = {};
    end

    vmcboundary.lightsource = ...
        extendCellArray(vmcboundary.lightsource, size(vmcmesh.BH,1));

    % build exterior refractive index
    if(exist('vmcmedium_in'))
       if(isfield(vmcmedium_in, 'refractive_index'))
           % avoid multidimensional indexing
           refractive_index = duplicateArray(vmcmedium_in.refractive_index(:), size(vmcmesh.H,1));
           % do not overwrite an existing exterior refractive index
           if(exist('vmcboundary_in') && ~isfield(vmcboundary_in, 'exterior_refractive_index'))
               s=findBoundaryElementsFromElements(vmcmesh.H, vmcmesh.BH);
               vmcboundary.exterior_refractive_index = refractive_index(s);
            end
       else
           warning(['Medium has no refractive index, hence exterior ' ...
                    'refrative index cannot be constructed']);
       end
    end

    if(exist('vmcboundary_in'))
        if(isfield(vmcboundary_in, 'lightsource'))
            vmcboundary.lightsource = ...
                extendCellArray(vmcboundary_in.lightsource, size(vmcmesh.BH,1));
        end
        if(isfield(vmcboundary_in, 'lightsource_direction'))
            vmcboundary.lightsource_direction = ...
                duplicateArray(vmcboundary_in.lightsource_direction, size(vmcmesh.BH,1));
        end

        if(isfield(vmcboundary_in, 'lightsource_position'))
            vmcboundary.lightsource_position = ...
                duplicateArray(vmcboundary_in.lightsource_position, size(vmcmesh.BH,1));
        end

        if(isfield(vmcboundary_in, 'exterior_refractive_index'))
            vmcboundary.exterior_refractive_index = ...
                duplicateArray(vmcboundary_in.exterior_refractive_index, size(vmcmesh.BH,1));
        end

        if(isfield(vmcboundary_in, 'lightsource_direction_type'))
            vmcboundary.lightsource_direction_type = ...
                   extendCellArray(vmcboundary_in.lightsource_direction_type, size(vmcmesh.BH,1));
        else
% force light source direction type to exist
             vmcboundary_in.lightsource_direction_type = {};
	         vmcboundary.lightsource_direction_type = ...
                extendCellArray(vmcboundary_in.lightsource_direction_type, size(vmcmesh.BH,1));
	end

        if(isfield(vmcboundary_in, 'lightsource_irradiance'))
            vmcboundary.lightsource_irradiance = ...
                duplicateArray(vmcboundary_in.lightsource_irradiance, size(vmcmesh.BH,1));
        end

        if(isfield(vmcboundary_in, 'lightsource_gaussian_sigma'))
            vmcboundary.lightsource_gaussian_sigma = ...
                duplicateArray(vmcboundary_in.lightsource_gaussian_sigma, size(vmcmesh.BH,1));
        end;

    end
end


function array_out = duplicateArray(array_in, desired_size, defaultvalue)
%
% function array_out = duplicateArray(array_in, size)
%
% INPUT
%
%  array_in:            incomplete array
%  desired_size:        desired size for the array
%
% OUTPUT
%
%  array_out:           completed array [size]
%
   % convert to column vector
   if(size(array_in,2) > size(array_in,1))
       array_in = transpose(array_in);
   end

   if(size(array_in,1) > desired_size)
      tmp = array_in(1:desired_size);
   else
      tmp = repmat(array_in,ceil(desired_size/size(array_in,1)),1);
   end

   array_out = tmp(1:desired_size,:);

end

function array_out  = extendCellArray(array_in, desired_size)
%COMPLETECELLARRAY Changes the size a cell array
%
% DESCRIPTION:
%
%       A helper function to extend a cell array
%
% USAGE:
%
% INPUTS:
%
%       array_in      Input array
%       desired_size  Desired size for the array
%
%
% OUTPUTS:
%
%       array_out      Extended array

  if(length(array_in) > desired_size)
     warning('extendCellArray is intended only for extending arrays.');
  end
  array_out = cell(desired_size, 1);
  array_out(1:length(array_in)) = array_in(1:length(array_in));
end


function output = findBoundaryElementsFromElements(H, BH)   
   sortedBH = transpose(sort(transpose(BH)));
   sortedH = transpose(sort(transpose(H)));
   if(size(H,2) == 4) 
       [test loc1] =ismember(sortedBH,sortedH(:,[1 2 3]),'rows'); 
       [test2 loc2] =ismember(sortedBH,sortedH(:,[1 2 4]),'rows'); 
       [test3 loc3] =ismember(sortedBH,sortedH(:,[1 3 4]),'rows'); 
       [test4 loc4] =ismember(sortedBH,sortedH(:,[2 3 4]),'rows'); 
       output=max([loc1 loc2 loc3 loc4]');
   else
       [test loc1] =ismember(sortedBH,sortedH(:,[1 2]),'rows'); 
       [test2 loc2] =ismember(sortedBH,sortedH(:,[1 3]),'rows'); 
       [test3 loc3] =ismember(sortedBH,sortedH(:,[2 3]),'rows'); 
       output=max([loc1 loc2 loc3]');
   end
end     

