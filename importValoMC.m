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
        timedep = fscanf(fp, '%e %e\n', 2);
        options.f = timedep(1);
        options.phase0 = timedep(2);
                
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
        if(dimension == 2) 
            vmcmesh.H = fscanf(fp, '%d %d %d\n', [3 Ne])' + 1;
            vmcmesh.BH = fscanf(fp, '%d %d\n', [2 Nb])' + 1;
            vmcmesh.r = fscanf(fp, '%e %e\n', [2 Nr])';
        else
            vmcmesh.H = fscanf(fp, '%d %d %d %d\n', [4 Ne])' + 1;
            vmcmesh.BH = fscanf(fp, '%d %d %d\n', [3 Nb])' + 1;
            vmcmesh.r = fscanf(fp, '%e %e %e\n', [3 Nr])';            
        end
       
        A = fscanf(fp, '%e %e %e %e\n', [4 Ne])';
      
        vmcmedium.absorption_coefficient = A(:, 1);
        vmcmedium.scattering_coefficient = A(:, 2);
        vmcmedium.scattering_anisotropy = A(:, 3);
        vmcmedium.refractive_index = A(:, 4);

        BCType = fscanf(fp, '%c\n', Nb);
        
        if(~feof(fp))
            BCn = fscanf(fp, '%e\n', Nb);
            if(~feof(fp))
                if(dimension == 2)
                   BCLNormal = fscanf(fp, '%e %e\n', 2 * Nb)';
                else
                   BCLNormal = fscanf(fp, '%e %e %e\n', 3 * Nb)';
                end
            end;
            if(~feof(fp))
                BCLightDirectionType = fscanf(fp, '%c\n', Nb)';
            end
            if(~feof(fp))
                BCLIntensity = fscanf(fp, '%e\n', Nb)';
            end
            if(~feof(fp))
                if(dimension == 2)
                   GaussianSigma = fscanf(fp, '%e\n', Nb)';
                end
            end
        end;
        fclose(fp);
    end
    
    vmcboundary.lightsource = characterToControlString(BCType,1);
    
    if(exist('BCn'))
       vmcboundary.exterior_refractive_index = BCn;    
    end

    if(exist('BCLNormal') && length(BCLNormal))
       vmcboundary.lightsource_direction = BCLNormal;    
    end
    
    if(exist('BCLightDirectionType')  && length(BCLightDirectionType))
       vmcboundary.lightsource_direction_type = characterToControlString(BCType,2);    
    end

    if(exist('GaussianSigma')  && length(GaussianSigma))
       vmcboundary.lightsource_gaussian_sigma = GaussianSigma;    
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
