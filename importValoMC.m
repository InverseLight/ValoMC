function [vmcmesh vmcmedium vmcboundary options solution] = ...
        importValoMC(problemdefinition_filename, solution_filename)

    dimension=getInputDimension(problemdefinition_filename);

    if(dimension ~= 2  && dimension ~= 3)
        error('Could not read input.');
    end

    if(~isempty(problemdefinition_filename));

        fp = fopen(problemdefinition_filename, 'r');
        A = fscanf(fp, '%d %d %d %d\n', 4);
        Ne = A(1); Nb = A(2); Nr = A(3);

        options.photon_count = A(4);
        timedep = fscanf(fp, '%e %e', 2);
        options.f = timedep(1);
        options.phase0 = timedep(2);
        timedep = fscanf(fp, '%d %d\n', 2);
        options.seed = timedep(1);
        solution.seed_used = options.seed;

        % medium dimensions for grid
        if(dimension == 2)
            B=fscanf(fp, '%i %i\n', 2);
            nx=B(1);
            ny=B(2);
            if(nx > 0 && ny > 0)
                vmcmedium.nx = nx;
                vmcmedium.ny = ny;
            end
        else
            B=fscanf(fp, '%i %i %i\n', 3);
            nx=B(1);
            ny=B(2);
            nz=B(3);
            if(nx > 0 && ny > 0 && nz > 0)
                vmcmedium.nx = nx;
                vmcmedium.ny = ny;
                vmcmedium.nz = nz;
            end
        end
    
        while(~feof(fp))
            tline = fgetl(fp);
            if(dimension == 2)
                if(strcmp(tline, 'H'))
                    vmcmesh.H = getArrayFromFile(fp, '%d %d %d\n', Ne)+1;  % not the same for 3d and 2d
                elseif(strcmp(tline, 'BH'))
                    vmcmesh.BH = getArrayFromFile(fp, '%d %d\n', Nb)+1;  % not the same for 3d and 2d
                elseif(strcmp(tline, 'r'))
                    vmcmesh.r = getArrayFromFile(fp, '%e %e\n', Nr);  % not the same for 3d and 2d
                elseif(strcmp(tline, 'BCLightDirection'))
                    vmcboundary.lightsource_direction = getArrayFromFile(fp, '%e %e\n', Nb);  % not the same for 3d and 2d  
                end                   
            else
                if(strcmp(tline, 'H'))
                    vmcmesh.H = getArrayFromFile(fp, '%d %d %d %d\n', Ne)+1;  % not the same for 3d and 2d
                elseif(strcmp(tline, 'BH'))
                    vmcmesh.BH = getArrayFromFile(fp, '%d %d %d\n', Nb)+1;  % not the same for 3d and 2d
                elseif(strcmp(tline, 'r'))
                    vmcmesh.r = getArrayFromFile(fp, '%e %e %e\n', Nr);  % not the same for 3d and 2d
                elseif(strcmp(tline, 'BCLightDirection'))
                    vmcboundary.lightsource_direction = getArrayFromFile(fp, '%e %e %e\n', Nb);  % not the same for 3d and 2d   
                end                  
            end
            if(strcmp(tline, 'mua mus g n'))
                A = getArrayFromFile(fp, '%e %e %e %e', Ne);
                vmcmedium.absorption_coefficient = A(:, 1);
                vmcmedium.scattering_coefficient = A(:, 2);
                vmcmedium.scattering_anisotropy = A(:, 3);
                vmcmedium.refractive_index = A(:, 4);
            elseif(strcmp(tline, 'BCType'))
                vmcboundary.lightsource = characterToControlString(getArrayFromFile(fp, '%c\n', Nb), 1); 
            elseif(strcmp(tline, 'BCn'))
                vmcboundary.exterior_refractive_index = getArrayFromFile(fp, '%e', Nb);                         
            elseif(strcmp(tline, 'BCLightDirectionType'))
                vmcboundary.lightsource_direction_type = characterToControlString(getArrayFromFile(fp, '%c', Nb), 2);
            elseif(strcmp(tline, 'BCLIntensity'))
                vmcboundary.lightsource_irradiance = getArrayFromFile(fp, '%e', Nb);
            elseif(strcmp(tline, 'GaussianSigma'))
                vmcboundary.lightsource_gaussian_sigma  = getArrayFromFile(fp, '%e', Nb);
            end
        end

        if(nargin == 2)
            if(~isempty(solution_filename))
                tmp = load(solution_filename);
                tmp = tmp(:, 1) + j * tmp(:, 2);
                solution.simulation_time = tmp(1);
                solution.element_fluence = tmp(2:size(vmcmesh.H,1)+1);
                solution.boundary_exitance = tmp(size(vmcmesh.H,1)+2:end);    
                if(dimension == 2)
                    if(isfield(vmcmedium,'nx') && isfield(vmcmedium,'ny'))
                        % Two dimensional input
                        first = reshape(solution.element_fluence(1:length(vmcmesh.H)/2),vmcmedium.nx, vmcmedium.ny);
                        second = reshape(solution.element_fluence(length(vmcmesh.H)/2+1:length(vmcmesh.H)),vmcmedium.nx, vmcmedium.ny);
                        solution.grid_fluence = (first+second)*0.5;
                    end
                else
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
        end
    end
end


function c = characterToControlString(str, mode)
%CHARACTERTOCONTROLSTRING Converts characters used in the C++ program into strings    
    c = {};
    if(mode == 1)
        % lightsources
        for k = 1:length(str)
            if(str(k) == 'a') c(k) = {'none'}; end;
            if(str(k) == 'l') c(k) = {'direct'}; end;
            if(str(k) == 'c') 
                c(k) = {'cosinic'}; 
            end;
            if(str(k) == 'g') c(k) = {'gaussian'}; end;
            if(str(k) == 'i') c(k) = {'isotropic'}; end;
            if(str(k) == 'p') c(k) = {'pencil'}; end;                        
        end
        if(mode == 2)
            % light source directions
            if(str(k) == 'a') c(k) = {'absolute'}; end;
            if(str(k) == 'r') c(k) = {'relative'}; end;
        end
    end
end

function array = getArrayFromFile(fp, formatspec, isize)
% GETARRAYFROMFILE Reads an array from an open file

% read the first line to obtain the dimensions
   tline = fgetl(fp);
   firstline = sscanf(tline, formatspec, isize)';
   jsize = size(firstline,2);
   array = zeros(isize, jsize);
   array(1, :) = firstline;
% now read the rest
   for ii=2:isize
      tline = fgetl(fp);
      array(ii, :) = sscanf(tline, formatspec, isize)';
   end

end
