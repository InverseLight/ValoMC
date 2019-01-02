function exportMMC(vmcmesh, mmc_suffix)
%EXPORTMMC Export the mesh for mesh-based Monte Carlo (MMC) software suite
%
% DESCRIPTION:
%
%       Exports the mesh for mesh-based Monte Carlo by Fang et
%       al.
%
% NOTE:
%
%       This is a simple wrapper but could be extended in the future to
%       export also the optical coefficients if there is a need for such
%       functionality
%
% USAGE:
%
%       exportMMC(vmcmesh, vmcmedium, mmc_suffix)
%       exportMMC(vmcmesh, vmcmedium, 'mymesh')
%
% INPUT:
%
%       vmcmesh        - mesh structure for ValoMC
%       vmcmedium      - medium for ValoMC
%       mmc_suffix     - prefix to be used in MMC, e.g. node_{prefix}.dat
%
% This function provided with ValoMC

   mmcinstalled = exist('savemmcmesh');
   if(~mmcinstalled)
       error('Could not detect MMC')
   end
   if(size(vmcmesh.H,2) ~= 4)
       error('2d not supported')
   end
   savemmcmesh(mmc_suffix,vmcmesh.r,vmcmesh.H);

end
