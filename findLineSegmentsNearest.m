function [segments] = findLineSegmentsNearest(vmcmesh, locations)
%
% Finds line segments nearest to each position
%
% INPUT
%
%  vmcmesh:         (described in documentation/list of structures)
%  locations:    an array that contains a position vector in each row.
%
% OUTPUT
%
%  segments:     the indices of the line segments nearest to the positions
%

    avgx = (vmcmesh.r(vmcmesh.BH(:,1),1) + vmcmesh.r(vmcmesh.BH(:,2),1))/2.0;
    avgy = (vmcmesh.r(vmcmesh.BH(:,1),2) + vmcmesh.r(vmcmesh.BH(:,2),2))/2.0;
    pos = [avgx avgy];
    segments = zeros(size(locations,1),1);
    for ii=1:size(locations,1)
       m=(pos - locations(ii,:)) .^2;
       norms=sum(m')';
       [minvalue minindex] = min(norms);
       segments(ii) = minindex; 
    end

end


