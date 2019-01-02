function vmcmesh = createCircularMesh(radius, dr)
%CREATECIRCULARMESH Creates a circular finite element mesh
%
% USAGE:
%
%       vmcmesh = createCircularMesh(radius, dh)
%
% DESCRIPTION:
%
%       Creates a circular finite element mesh (2D)
%
% INPUT:
%
%       radius      - Total radius of the circle [mm]
%       dr          - Separation of each node in radius [mm]
%
% OUTPUT:
%
%       vmcmesh     - Resulting circular mesh
%
% This function is provided with ValoMC

    RAD = radius;
    vmcmesh.r = [];
    for r = linspace(0, RAD, ceil(RAD / dr))
        Nth = max([1, ceil(2 * pi * r / dr)]);
        th = 2 * pi * (0 : Nth - 1)' / Nth;
        vmcmesh.r = [ vmcmesh.r ; r * cos(th), r * sin(th) ];
    end;

    %bndr = boundary(vmcmesh.r);
    %vmcmesh.BH = ones(length(bndr)-1,2);

    %for i = 1 : length(bndr)-1
    %   vmcmesh.BH(i,1)=bndr(i);
    %   vmcmesh.BH(i,2)=bndr(i+1);
    %end

    vmcmesh.H = delaunay(vmcmesh.r(:, 1), vmcmesh.r(:, 2));
    vmcmesh.BH = createBH(vmcmesh.H, createHN(vmcmesh.H));

end

