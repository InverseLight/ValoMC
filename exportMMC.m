function exportMMC(vmcmesh, vmcmedium, mmc_suffix)
%EXPORTMMC Export the simulation setup for mesh-based Monte Carlo
%
% DESCRIPTION:
%
%       Exports the simulation setup for mesh-based Monte Carlo by Fang et
%       al. Note that only the medium properties are set.
%       
% USAGE:
%
%       exportX3D(vmcmesh, vmcmedium, mmc_suffix)
%       exportX3D(vmcmesh, vmcmedium, 'mymesh')
%
% INPUTS:
%
%       vmcmesh        - mesh structure for ValoMC
%       vmcmedium      - medium for ValoMC
%       mmc_suffix     - prefix to be used in MMC, e.g. node_{prefix}.dat
%
   mmcinstalled = exist('savemmcmesh');
   if(~mmcinstalled)
       error('Could not detect MMC')
   end
   if(size(vmcmesh.H,2) ~= 4)
       error('2d not supported')
   end
   savemmcmesh(mmc_suffix,vmcmesh.r,vmcmesh.H);
   
end
