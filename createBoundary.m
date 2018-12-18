function [vmcboundary, vmcmesh_out] = createBoundary(vmcmesh, vmcmedium_in, vmcboundary_in)
% CREATEBOUNDARY Create a boundary structure for a mesh
%
% DESCRIPTION:
%
%       This function is used to create a boundary structure.
%       The functionality depends on the number of input arguments
%
%       1: vmcmesh: create a minimal number of fields in order to run
%       a simulation, i.e. a lightsource field whose length is
%       equal to the number of boundary elements
%
%       2: vmcmesh, vmcmedium: in addition to 1, create
%       exterior_refactive_index so that photons experience no change
%       in refractive index as they leave the medium
%
%       3: vmcmesh, vmcmedium, vmcboundary: repeat any entries in
%       vmcboundary so that the lengths of the arrays are equal to the
%       number of boundary elements in vmcmesh. If no
%       exterior_refractive_index exists in vmcboundary, creates one
%       according section 2.
%
%  USAGE:
%        vmcboundary = createBoundary(vmcmesh)
%        vmcboundary = createBoundary(vmcmesh, vmcmedium)
%        vmcboundary = createBoundary(vmcmesh, vmcmedium, vmcboundary)
%        [vmcboundary, vmcmesh] = createBoundary(vmcmesh)
%
% INPUTS:
%       vmcmesh      - mesh structure
%
% OPTIONAL INPUTS:
%
%       vmcmedium    - medium structure
%       vmcboundary  - vmcboundary
%
% OUTPUT:
%
%       vmcboundary - boundary constructed in the way explained in description
%
% OPTIONAL OUTPUT:
%
% vmcmesh:
%       vmcmesh     - if vmcmesh.BH is missing from the input, a structure
%                     that has it
%
    if(nargout == 2)
        vmcmesh_out = vmcmesh;
        if(~isfield(vmcmesh, 'HN'))
           vmcmesh_out.HN = createBH(vmcmesh.H,createHN(vmcmesh.H));
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
           if(~isfield(vmcmedium_in, 'exterior_refractive_index'))
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

