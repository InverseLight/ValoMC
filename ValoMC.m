function solution = ValoMC(vmcmesh, vmcmedium, vmcboundary, vmcoptions)
% VALOMC Runs a photon transport simulation
%
% USAGE:
%
%       solution = ValoMC(vmcmesh, vmcmedium, vmcboundary)
%       solution = ValoMC(vmcmesh, vmcmedium, vmcboundary, vmcoptions)
%
% DESCRIPTION:
%
%       This main function used to initiate simulations.
%       The input structures are documented in more detail
%       in the homepage (see below).
%
% INPUT:
%
%       vmcmesh                   - contains the geometry of the system
%         .H                      - element topology
%         .BH                     - boundary element topology
%         .r                      - node coordinates
%       vmcmedium                 - contains the optical properties of the medium
%         .absorption_coefficient   all have same size as size(H,1)
%         .scattering_coefficient
%         .scattering_anisotropy
%         .refractive_index
%       vmcboundary
%         .lightsource            - type of the lightsource (e.g. {'cosinic'}. {'gaussian'})
%         .lightsource_direction  - size(H, Ndim)
%         .lightsource_direction_type  - {'relative'} or {'absolute'}, direction with respect to normal or an absolute direction
%
% OPTIONAL INPUT:
%       vmcoptions
%         .photon_count
%         .disable_progressbar    - true or false
%         .seed                   - random number generator seed
%
%
% OUTPUT:
%
%       output
%         .element_fluence         - size(H,1)
%         .boundary_exitance       - size(BH,1)
%         .seed_used
%         .simulation_time
%
% OPTIONAL OUTPUT
%
%         .grid_fluece             - per pixel/voxel fluence if the mesh was created using createGridMesh
%
% This function is provided with ValoMC


% VALOMC Runs a photon transport simulation.
%
% INPUT
%
%   vmcmesh   : 
%
% OUTPUT
%
%          see documentation/structure reference
%
% SEE ALSO
%
%





% Gather the input for the mex code.
%
% Structure of input:
%          H   Triangular topology matrix [Ne x 3]
%         HN   Topology neighbourhood (optional) [Ne x 3]
%         BH   Boundary topology of H [Nb x 2]
%    BDomain   Subdomains for elements in BH [Nb]
%          r   Grid nodes (mm) [Np x 2]
%        mua   Absorption coefficient (1/mm) for each Domain
%        mus   Scattering coefficient (1/mm) for each Domain
%          g   Scattering anisotropy for each Domain
%          n   Index of refraction for each Domain
%     BCType   Boundary condions for each BDomain
%  BCLNormal   Direction of light source for each BDomain (optional)
% GaussianSigma   Parameter for Gaussian light source (optional)
%        BCn   Exterior index of refraction for each BDomain (optional)
%          f   Frequency of modulation for the light source
%    Nphoton   Number of photons to compute
%     rnseed   Random number generator seed

    if (nargin < 3)
        error('Not enought input arguments.');
    end

    if(~isfield(vmcmesh, 'H'))
        error('Mesh does not contain element topology')
    end

    dimensionality = 0;

    if(size(vmcmesh.H,2) == 3)
        dimensionality = 2;
    else
        dimensionality = 3;
    end
    
    
    % Add optional fields related to vmcmedium
    
    vmcmedium = createMedium(vmcmesh,vmcmedium);
    
    % Add optional fields related to boundary

    vmcboundary = createBoundary(vmcmesh, vmcmedium, vmcboundary);
    
    % Construct exterior index of refraction but avoid overriding
    if(~isfield(vmcboundary, 'exterior_refractive_index'))
        vmcboundary_ext = createBoundary(vmcmesh, vmcmedium);
        vmcboundary.exterior_refractive_index = vmcboundary_ext.exterior_refractive_index;
    end
    
    if(~isfield(vmcmesh, 'HN'))
        vmcmesh.HN = [];
    end

    if(isfield(vmcmedium,'ny') || isfield(vmcmedium,'nz'))
        vmcmedium.absorption_coefficient = duplicateArray(vmcmedium.absorption_coefficient(:), length(vmcmesh.H));
    end
    if(isfield(vmcmedium,'ny') || isfield(vmcmedium,'nz'))
        vmcmedium.scattering_coefficient = duplicateArray(vmcmedium.scattering_coefficient(:), length(vmcmesh.H));
    end
    if(isfield(vmcmedium,'ny') || isfield(vmcmedium,'nz'))
        vmcmedium.scattering_anisotropy = duplicateArray(vmcmedium.scattering_anisotropy(:), length(vmcmesh.H));
    end
    if(isfield(vmcmedium,'ny') || isfield(vmcmedium,'nz'))
        vmcmedium.refractive_index = duplicateArray(vmcmedium.refractive_index(:), length(vmcmesh.H));
    end

    % Convert fields
    H = int64(vmcmesh.H - 1);
    HN = int64(vmcmesh.HN);
    BH = int64(vmcmesh.BH - 1);
    r = double(vmcmesh.r);
    mua = double(vmcmedium.absorption_coefficient);
    mus = double(vmcmedium.scattering_coefficient);
    g = double(vmcmedium.scattering_anisotropy);
    n = double(vmcmedium.refractive_index);
    BCType = int8(arrayfun(@(x) controlStringToCharacter(x,'a'), vmcboundary.lightsource));

    % Handle optional fields
    if(~isfield(vmcboundary,'lightsource_irradiance'))
        vmcboundary.lightsource_irradiance = [];
    end
    BCIntensity = double(vmcboundary.lightsource_irradiance);
    if(~isfield(vmcboundary,'lightsource_direction'))
        vmcboundary.lightsource_direction = [];
    end
    if(isfield(vmcboundary,'lightsource_gaussian_sigma'))
        GaussianSigma = double(vmcboundary.lightsource_gaussian_sigma);
    end
    if(~isfield(vmcboundary,'lightsource_direction_type'))
        error('Lightsource direction type is not set.');
    else
        % make empty BCLightDirectionTypes 'n' so that it won't be used
        BCLightDirectionType = int8(arrayfun(@(x) controlStringToCharacter(x,'n'), vmcboundary.lightsource_direction_type));
    end

    BCLightDirection = double(vmcboundary.lightsource_direction);
    BCn = double(vmcboundary.exterior_refractive_index);

    % set default options
    f = double(0.0);
    Nphoton = int64(1e6);
    phase0 = 0;
    disable_pbar = int64(0);

    % complement with user provided options
    if(exist('vmcoptions')==1)
        if(isfield(vmcoptions, 'frequency'))
            f = double(vmcoptions.frequency);
        end
        if(isfield(vmcoptions,'photon_count'))
            Nphoton = int64(vmcoptions.photon_count);
        end
        if(isfield(vmcoptions, 'phase0'))
            phase0 = double(vmcoptions.phase0);
        end
        if(isfield(vmcoptions, 'disable_progressbar'))
            if(vmcoptions.disable_progressbar) 
                disable_pbar = int64(1);
            end
        end
	if(isfield(vmcoptions, 'seed'))
	    rnseed(1) = vmcoptions.seed;
	    rnseed(2) = 1;
	end
    else
        vmcoptions = struct();
    end

    if(~exist('rnseed'))
       rnseed(1) = 0;
       rnseed(2) = 0;
    end

    if(dimensionality == 2)
        if(any(find(BCType == int8('p'))))
            error('Pencil beam currently not supported in 2D.');
        end


        if(isfield(vmcoptions,'export_filename'))
            fp = fopen(vmcoptions.export_filename, 'w');
            fprintf(fp, '%d %d %d %d %d %d\n', size(H, 1), size(BH, 1), size(r, ...
                                                              1), Nphoton);
            fprintf(fp, '%e %e %d %d\n', f, phase0, rnseed(1), rnseed(2));

            % write the array dimensions to the file so that import
            % can know them

            if(isfield(vmcmedium,'nx') && isfield(vmcmedium,'ny'))
                fprintf(fp, '%i %i\n',vmcmedium.nx, vmcmedium.ny);
            else
                fprintf(fp, '0 0\n');
            end
            fprintf(fp, 'H\n');
            fprintf(fp, '%d %d %d\n', H');
            fprintf(fp, 'BH\n');
            fprintf(fp, '%d %d\n', BH');
            fprintf(fp, 'r\n');
            fprintf(fp, '%e %e\n', r');
            fprintf(fp, 'mua mus g n\n');
            fprintf(fp, '%e %e %e %e\n', [ mua(:) mus(:) g(:) n(:) ]');
            fprintf(fp, 'BCType\n');
            fprintf(fp, '%c\n', BCType);

            if exist('BCn')
                fprintf(fp, 'BCn\n');
                fprintf(fp, '%e\n', BCn);
            end
            if exist('BCLightDirection')
                fprintf(fp, 'BCLightDirection\n');
                fprintf(fp, '%e %e\n', BCLightDirection');
            end
            if exist('BCLightDirectionType')
                fprintf(fp, 'BCLightDirectionType\n');
                fprintf(fp, '%c\n', BCLightDirectionType');
            end
            if exist('BCLIntensity')
                fprintf(fp, 'BCLIntesity\n');
                fprintf(fp, '%e\n', BCLIntensity');
            end
            if exist('GaussianSigma')
                fprintf(fp, 'GaussianSigma\n');
                fprintf(fp, '%e\n', GaussianSigma');
            end

            fclose(fp);

            return
        else

        if(~exist('GaussianSigma'))
           GaussianSigma = [];
        end
        % Solve
        [solution.element_fluence, solution.boundary_exitance, solution.boundary_fluence, solution.simulation_time, solution.seed_used] = MC2Dmex(H, HN, BH, r, BCType, BCIntensity, BCLightDirectionType, BCLightDirection, BCn, mua, mus, g, n, f, phase0, Nphoton, GaussianSigma, disable_pbar, int64(rnseed));
            if(isfield(vmcmedium,'nx') && isfield(vmcmedium,'ny'))
                % Two dimensional input
                first = reshape(solution.element_fluence(1:length(vmcmesh.H)/2),vmcmedium.nx, vmcmedium.ny);
                second = reshape(solution.element_fluence(length(vmcmesh.H)/2+1:length(vmcmesh.H)),vmcmedium.nx, vmcmedium.ny);
                solution.grid_fluence = (first+second)*0.5;
            end
        end
      else if(dimensionality == 3)
        if(any(find(BCType == int8('p'))))
            if(~isfield(vmcboundary,'lightsource_position') || length(vmcboundary.lightsource_position) ~= size(BH,1)) 
                error('Please provide relative positions for pencil beam using lightsource_position');
            end
            vmcboundary.lightsource_direction_type = ...
                extendCellArray(vmcboundary.lightsource_direction_type, size(vmcmesh.BH,1));
            BCLightDirectionType = int8(arrayfun(@(x) controlStringToCharacter(x,'n'), vmcboundary.lightsource_direction_type));

            % if no light direction was provided previously, overwrite it with position
            if(size(BCLightDirection,1) ~= size(BH,1))
                BCLightDirection = double(vmcboundary.lightsource_position);
            else
                % overwrite only elements with a pencil lightsource
                BCLightDirection(find(BCType == int8('p')),1) = double(vmcboundary.lightsource_position(find(BCType == int8('p'))), 1);
                BCLightDirection(find(BCType == int8('p')),2) = double(vmcboundary.lightsource_position(find(BCType == int8('p'))), 2);
                BCLightDirection(find(BCType == int8('p')),3) = double(vmcboundary.lightsource_position(find(BCType == int8('p'))), 3);
            end
            BCLightDirectionType(find(BCType == int8('p'))) = int8('p');
        end
        if(isfield(vmcoptions,'export_filename'))        
            
            fp = fopen(vmcoptions.export_filename, 'w');
            fprintf(fp, '%d %d %d %d\n', size(H, 1), size(BH, 1), size(r, 1), Nphoton);    
            fprintf(fp, '%18.10f %18.10f %d %d\n', f, phase0, rnseed(1), rnseed(2));

            if(isfield(vmcmedium,'nx') && isfield(vmcmedium,'ny') && ...
               isfield(vmcmedium,'nz'))
                fprintf(fp, '%i %i %i\n', vmcmedium.nx, vmcmedium.ny, vmcmedium.nz);
            else 
                fprintf(fp, '0 0 0\n');
            end 
            fprintf(fp, 'H\n');
            fprintf(fp, '%d %d %d %d\n', H');
            fprintf(fp, 'BH\n');    
            fprintf(fp, '%d %d %d\n', BH');
            fprintf(fp, 'r\n');
            fprintf(fp, '%18.10f %18.10f %18.10f\n', r');
            fprintf(fp, 'mua mus g n\n');
            fprintf(fp, '%18.10f %18.10f %18.10f %18.10f\n', [ mua(:) mus(:) g(:) n(:) ]');
            fprintf(fp, 'BCType\n');
            fprintf(fp, '%c\n', BCType);
            if exist('BCn')
                fprintf(fp, 'BCn\n');
                fprintf(fp, '%18.10f\n', BCn);
            end
            if exist('BCLightDirection')
                fprintf(fp, 'BCLightDirection\n');
                fprintf(fp, '%18.10f %18.10f %18.10f\n', BCLightDirection');
            end
            if exist('BCLightDirectionType')
                fprintf(fp, 'BCLightDirectionType\n');
                fprintf(fp, '%c\n', char(BCLightDirectionType'));
            end
            if exist('BCLIntensity')
                fprintf(fp, 'BCLIntesity\n');
                fprintf(fp, '%18.10f\n', BCLIntensity');
            end
            if exist('GaussianSigma')
                fprintf(fp, 'GaussianSigma\n');
                fprintf(fp, '%18.10f\n', GaussianSigma');
            end
            fclose(fp);
            return
        else
            [solution.element_fluence, solution.boundary_exitance, solution.boundary_fluence, solution.simulation_time, solution.seed_used] = MC3Dmex(H, HN, BH, r, BCType, BCIntensity, BCLightDirectionType, BCLightDirection, BCn, mua, mus, g, n, f, phase0, Nphoton,disable_pbar, int64(rnseed));
        end
        if(isfield(vmcmedium,'nx') && isfield(vmcmedium,'ny') && isfield(vmcmedium,'nz'))
            % Three dimensional input
            nvoxels= length(vmcmesh.H) / 6;
            first = reshape(solution.element_fluence(1:nvoxels), vmcmedium.nx, vmcmedium.ny, vmcmedium.nz);
            second = reshape(solution.element_fluence(nvoxels+1:2*nvoxels),vmcmedium.nx, vmcmedium.ny, vmcmedium.nz);
            third = reshape(solution.element_fluence(nvoxels*2+1:3*nvoxels),vmcmedium.nx, vmcmedium.ny, vmcmedium.nz);
            fourth = reshape(solution.element_fluence(nvoxels*3+1:4*nvoxels),vmcmedium.nx, vmcmedium.ny, vmcmedium.nz);
            fifth = reshape(solution.element_fluence(nvoxels*4+1:5*nvoxels),vmcmedium.nx, vmcmedium.ny, vmcmedium.nz);
            sixth = reshape(solution.element_fluence(nvoxels*5+1:6*nvoxels),vmcmedium.nx, vmcmedium.ny, vmcmedium.nz);

            solution.grid_fluence = (first+second+third+fourth+fifth+sixth)/6;
        end
      end
    end

    % Remove an imaginary solution that is zero
    if(f == 0)
        solution.element_fluence = real(solution.element_fluence);
        solution.boundary_fluence = real(solution.boundary_fluence);
        solution.boundary_exitance = real(solution.boundary_exitance);
    end

end

function array_out = duplicateArray(array_in, desired_size, defaultvalue)
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
%
% DESCRIPTION:
%       A helper function to quicky take take content of one
%       array and make another, bigger cell array with it.
%
% INPUTS:
%       array_in      Input array
%       desired_size  Desired size for the array
%
%
% OUTPUTS:
%       array_out      Extended array

  if(length(array_in) > desired_size)
     warning('extendCellArray is intended only for extending arrays.');
  end
  array_out = cell(desired_size, 1);
  array_out(1:length(array_in)) = array_in(1:length(array_in));
end

function c = controlStringToCharacter(string, default)
%CONTROLSTRINGTOCHARACTER Convert a control string to a single character for the C++ code 
%
%
% DESCRIPTION:
%       (ValoMC internal use only) 
%       Converts strings that are used in the ValoMC Matlab interface like 'cosinic'
%       into a single byte character. For example 'gaussian' -> 'g'
%   
% USAGE:
%       output = controlStringToCharacter('gaussian', 'n')
%
% INPUTS:
%       string     - Control string to convert
%       default    - If no character is found for conversion, return this value
%
% OUTPUTS:
%       c          - Converted character
    
      if(strcmp(string,'none')) c = 'a';
      elseif(strcmp(string,'direct')) c = 'l';
      elseif(strcmp(string,'cosinic')) c = 'c';
      elseif(strcmp(string,'gaussian')) c = 'g';
      elseif(strcmp(string,'isotropic')) c = 'i';
      elseif(strcmp(string,'absolute')) c = 'a';
      elseif(strcmp(string,'relative')) c = 'r';
      elseif(strcmp(string,'pencil')) c = 'p';
      else c = default;
      end
    
end
    
    


