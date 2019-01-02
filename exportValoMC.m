function exportValoMC(filename, vmcmesh, vmcmedium, vmcboundary, vmcoptions)
%EXPORTVALOMC Exports the simulation setup into ASCII file that can be run by the external executable
%
% USAGE:
%
%        exportValoMC(filename, vmcmesh, vmcmedium, vmcboundary, vmcoptions)
%
% DESCRIPTION:
%
%       Used to export the simulation setup into a file. The file can be opened by the external executable
%       and enables running the simulations without MATLAB e. g. on a computing cluster.
%
% INPUT:
%
%       vmcmesh                       - mesh structure, contains the geometry of the system
%       vmcmedium                     - medium structure, used to set e. g. optical coefficients
%       vmcboundary                   - boundary structure, used to set e. g. light sources
%       vmcoptions                    - contains additional options, e. g. number of photon packets
%
% EXAMPLE:
%
%       exportValoMC('problem_in.txt',vmcmesh, vmcmedium, vmcboundary, vmcoptions);
%       ... exit matlab run the simulation ....
%       MC2D problem_in.txt problem_out.txt
%       .. or ..
%       MC3D problem_in.txt problem_out.txt
%       ... return to matlab ...
%       [vmcmesh vmcmedium vmcboundary options solution] = importValoMC('problem_in.txt', 'problem_out.txt')
%
%
% This function is provided with ValoMC
    vmcoptions.export_filename = filename;
    ValoMC(vmcmesh, vmcmedium, vmcboundary, vmcoptions);
end
