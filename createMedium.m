function vmcmedium_out = createMedium(vmcmesh, vmcmedium)
%CREATEMEDIUM Creates a 'vmcmedium' structure array
%
% DESCRIPTION:
%        
%       Creates a 'vmcmedium' structure array (see https://inverselight.github.io/ValoMC/structures.html)
%       or turns the existing fields in a given medium into arrays that are compatible in size with a
%       given mesh.
%    
% USAGE:
%       vmcmedium = createMedium(vmcmesh);
%
%           returns the following structure
%    
%                   vmcmedium.refractive_index(:) = 1;               
%                   vmcmedium.scattering_coefficient(:) = 0; 
%                   vmcmedium.absorption_coefficient(:) = 0; 
%                   vmcmedium.scattering_anisotropy(:) = 0;
%
%           The size of each array is equal to the number of elements
%           in the mesh.
%
%       vmcmedium_out = createMedium(vmcmesh, vmcmedium);
%
%           returns the structure as above but each existing field in
%           vmcmedium is copied to vmcmedium_out by repeating it until
%           the size matches the number of elements in the mesh
%
% INPUTS:
%       vmcmesh       -  see https://inverselight.github.io/ValoMC/structures.html  
%
% OPTIONAL INPUTS:
%       vmcmedium     -  incomplete structure array for medium (see https://inverselight.github.io/ValoMC/structures.html)
%
% OUTPUTS:
%       vmcmedium_out -  structure array for medium

    if(~exist('vmcmedium'))
         vmcmedium.refractive_index = 1;               
         vmcmedium.scattering_coefficient = 0;
         vmcmedium.absorption_coefficient = 0; 
         vmcmedium.scattering_anisotropy = 0;        
    end
    vmcmedium_out = vmcmedium;

    if(~isfield(vmcmedium, 'refractive_index'))
       vmcmedium.refractive_index = 1;
    end;
    if(~isfield(vmcmedium, 'scattering_coefficient'))
         vmcmedium.scattering_coefficient = 0.0; 
    end;
    if(~isfield(vmcmedium, 'absorption_coefficient'))
         vmcmedium.absorption_coefficient = 0.0;
    end;
    if(~isfield(vmcmedium, 'scattering_anisotropy'))
         vmcmedium.absorption_coefficient = 0.0; 
    end;
    
    % check if the optical coefficients are given as multidimensional arrays

    if(~isvector(vmcmedium.refractive_index)) 
        vmcmedium_out.nx = size(vmcmedium.refractive_index,1);
        vmcmedium_out.ny = size(vmcmedium.refractive_index,2);
        if(size(vmcmedium.refractive_index, 3)>1)
		vmcmedium_out.nz = size(vmcmedium.refractive_index,3);        
           if(vmcmedium_out.nx*vmcmedium_out.ny*vmcmedium_out.nz ~= length(vmcmesh.H)/6)
               warning('Mesh size and medium size do not seem compatible.');
           end
        else
           if(vmcmedium_out.nx*vmcmedium_out.ny ~= length(vmcmesh.H)/2)
               warning('Mesh size and medium size do not seem compatible.');
           end
        end
    end
    
    if(~isvector(vmcmedium.scattering_coefficient))
        vmcmedium_out.nx = size(vmcmedium.scattering_coefficient,1);
        vmcmedium_out.ny = size(vmcmedium.scattering_coefficient,2);
        
        if(size(vmcmedium.scattering_coefficient, 3)>1)
		vmcmedium_out.nz = size(vmcmedium.scattering_coefficient,3);        
           if(vmcmedium_out.nx*vmcmedium_out.ny*vmcmedium_out.nz ~= length(vmcmesh.H)/6)
               warning('Mesh size and medium size do not seem compatible.');
           end
        else
           if(vmcmedium_out.nx*vmcmedium_out.ny ~= length(vmcmesh.H)/2)
               warning('Mesh size and medium size do not seem compatible.');
           end
        end
    end
    
    if(~isvector(vmcmedium.absorption_coefficient))
        vmcmedium_out.nx = size(vmcmedium.absorption_coefficient,1);
        vmcmedium_out.ny = size(vmcmedium.absorption_coefficient,2);
        if(size(vmcmedium.absorption_coefficient, 3)>1)
           vmcmedium_out.nz = size(vmcmedium.absorption_coefficient,3);
           if(vmcmedium_out.nx*vmcmedium_out.ny*vmcmedium_out.nz ~= length(vmcmesh.H)/6)
               warning('Mesh size and medium size do not seem compatible.');
           end
        else
           if(vmcmedium_out.nx*vmcmedium_out.ny ~= length(vmcmesh.H)/2)
               warning('Mesh size and medium size do not seem compatible.');
           end
        end
    end
    
    if(~isvector(vmcmedium.scattering_anisotropy))
        vmcmedium_out.nx = size(vmcmedium.scattering_anisotropy,1);
        vmcmedium_out.ny = size(vmcmedium.scattering_anisotropy,2);
        if(size(vmcmedium.scattering_anisotropy, 3)>1)
           vmcmedium_out.nz = size(vmcmedium.scattering_anisotropy,3);
           if(vmcmedium_out.nx*vmcmedium_out.ny*vmcmedium_out.nz ~= length(vmcmesh.H)/6)
               warning('Mesh size and medium size do not seem compatible.');
           end
        else
           if(vmcmedium_out.nx*vmcmedium_out.ny ~= length(vmcmesh.H)/2)
               warning('Mesh size and medium size do not seem compatible.');
           end
        end
    end


    vmcmedium_out.refractive_index = duplicateArray(vmcmedium.refractive_index(:), size(vmcmesh.H,1));
    vmcmedium_out.scattering_coefficient = duplicateArray(vmcmedium.scattering_coefficient(:), size(vmcmesh.H,1));
    vmcmedium_out.absorption_coefficient = duplicateArray(vmcmedium.absorption_coefficient(:), size(vmcmesh.H,1));
    vmcmedium_out.scattering_anisotropy = duplicateArray(vmcmedium.scattering_anisotropy(:), size(vmcmesh.H,1));
   
end
