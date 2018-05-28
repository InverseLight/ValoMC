function vmcmedium_out = createMedium(vmcmesh, vmcmedium_in)
%
% Complete a 'vmcmedium' -structure for a given geometry
%
% vmcmedium_out = createMedium(vmcmesh, vmcmedium_in)
%
%
% INPUT
%
%  vmcmesh:  	    mesh structure (see documentation/structure reference)
%
% OPTIONAL INPUT
%
%  vmcmedium_in:	existing medium structure (see documentation/structure reference)
%               where the array size does not match the number of elements in
%               the mesh.
% OUTPUT
%
%  vmcmedium_out:  vmcmedium structure for a given mesh. The array sizes match
%               the number of elements in the mesh.
%               If 'vmcmedium_in' was given as a input, the values for the arrays in
%               'vmcmedium_out' are obtained by repeating 'vmcmedium_in'. 
%               If no 'vmcmedium_in' was given, the function returns a structure with
%               the fields formatted as follows:
%
%                   vmcmedium_out(:).refractive_index = 1;               
%                   vmcmedium_out(:).scattering_coefficient = 0; [1/mm]
%                   vmcmedium_out(:).absorption_coefficient = 0; [1/mm]
%                   mediun_out(:).scattering_anisotropy = 0;
%
%
    if(~exist('vmcmedium_in'))
         vmcmedium_in.refractive_index = 1;               
         vmcmedium_in.scattering_coefficient = 0;
         vmcmedium_in.absorption_coefficient = 0; 
         vmcmedium_in.scattering_anisotropy = 0;        
    end
    vmcmedium_out = vmcmedium_in;

    if(~isfield(vmcmedium_in, 'refractive_index'))
       vmcmedium_in.refractive_index = 1;
    end;
    if(~isfield(vmcmedium_in, 'scattering_coefficient'))
         vmcmedium_in.scattering_coefficient = 0.0; 
    end;
    if(~isfield(vmcmedium_in, 'absorption_coefficient'))
         vmcmedium_in.absorption_coefficient = 0.0;
    end;
    if(~isfield(vmcmedium_in, 'scattering_anisotropy'))
         vmcmedium_in.absorption_coefficient = 0.0; 
    end;
    
    if(~isvector(vmcmedium_in.refractive_index)) 
        vmcmedium_out.nx = size(vmcmedium_in.refractive_index,1);
        vmcmedium_out.ny = size(vmcmedium_in.refractive_index,2);
        if(size(vmcmedium_in.refractive_index, 3)>1)
		vmcmedium_out.nz = size(vmcmedium_in.refractive_index,3);        
           if(vmcmedium_out.nx*vmcmedium_out.ny*vmcmedium_out.nz ~= length(vmcmesh.H)/6)
               warning('Mesh size and medium size do not seem compatible.');
           end
        else
           if(vmcmedium_out.nx*vmcmedium_out.ny ~= length(vmcmesh.H)/2)
               warning('Mesh size and medium size do not seem compatible.');
           end
        end
    end
    
    if(~isvector(vmcmedium_in.scattering_coefficient))
        vmcmedium_out.nx = size(vmcmedium_in.scattering_coefficient,1);
        vmcmedium_out.ny = size(vmcmedium_in.scattering_coefficient,2);
        
        if(size(vmcmedium_in.scattering_coefficient, 3)>1)
		vmcmedium_out.nz = size(vmcmedium_in.scattering_coefficient,3);        
           if(vmcmedium_out.nx*vmcmedium_out.ny*vmcmedium_out.nz ~= length(vmcmesh.H)/6)
               warning('Mesh size and medium size do not seem compatible.');
           end
        else
           if(vmcmedium_out.nx*vmcmedium_out.ny ~= length(vmcmesh.H)/2)
               warning('Mesh size and medium size do not seem compatible.');
           end
        end
    end
    
    if(~isvector(vmcmedium_in.absorption_coefficient))
        vmcmedium_out.nx = size(vmcmedium_in.absorption_coefficient,1);
        vmcmedium_out.ny = size(vmcmedium_in.absorption_coefficient,2);
        if(size(vmcmedium_in.absorption_coefficient, 3)>1)
           vmcmedium_out.nz = size(vmcmedium_in.absorption_coefficient,3);
           if(vmcmedium_out.nx*vmcmedium_out.ny*vmcmedium_out.nz ~= length(vmcmesh.H)/6)
               warning('Mesh size and medium size do not seem compatible.');
           end
        else
           if(vmcmedium_out.nx*vmcmedium_out.ny ~= length(vmcmesh.H)/2)
               warning('Mesh size and medium size do not seem compatible.');
           end
        end
    end
    
    if(~isvector(vmcmedium_in.scattering_anisotropy))
        vmcmedium_out.nx = size(vmcmedium_in.scattering_anisotropy,1);
        vmcmedium_out.ny = size(vmcmedium_in.scattering_anisotropy,2);
        if(size(vmcmedium_in.scattering_anisotropy, 3)>1)
           vmcmedium_out.nz = size(vmcmedium_in.scattering_anisotropy,3);
           if(vmcmedium_out.nx*vmcmedium_out.ny*vmcmedium_out.nz ~= length(vmcmesh.H)/6)
               warning('Mesh size and medium size do not seem compatible.');
           end
        else
           if(vmcmedium_out.nx*vmcmedium_out.ny ~= length(vmcmesh.H)/2)
               warning('Mesh size and medium size do not seem compatible.');
           end
        end
    end


    vmcmedium_out.refractive_index = duplicateArray(vmcmedium_in.refractive_index(:), size(vmcmesh.H,1));
    vmcmedium_out.scattering_coefficient = duplicateArray(vmcmedium_in.scattering_coefficient(:), size(vmcmesh.H,1));
    vmcmedium_out.absorption_coefficient = duplicateArray(vmcmedium_in.absorption_coefficient(:), size(vmcmesh.H,1));
    vmcmedium_out.scattering_anisotropy = duplicateArray(vmcmedium_in.scattering_anisotropy(:), size(vmcmesh.H,1));
   
end
