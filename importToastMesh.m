function vmcmesh = importToastMesh(toastMesh)
%IMPORTTOASTMESH Imports a mesh from the Toast++ software suite
%
% USAGE:
%
%       vmcmesh = importToastMesh(toastmesh)
%
% EXAMPLE:
%
%       rad = 25;
%       nsect = 6;
%       nring = 32;
%       nbnd = 2;
%       [vtx,idx,eltp] = mkcircle(rad,nsect,nring,nbnd);
%       toastmesh = toastMesh(vtx,idx,eltp);
%       vmcmesh = importToastMesh(toastmesh);
%
% DESCRIPTION:
%
%       This is a wrapper function that converts the data in
%       an instance of the toastMesh class into the mesh structure
%       of ValoMC.
%
% INPUT:
%
%       toastmesh   - toastMesh class
%
% OUTPUT:
%
%       vmcmesh     - vmcmesh structure
%
% This function is provided with ValoMC
%
    if not(isa(toastMesh, 'toastMesh'))
        error('Incorrect argument');
    end
    element=toastMesh.Element(1);
    if(element.elid ~= 1)
        error('Element type currently not supported');
    end
    [vmcmesh.r, vmcmesh.H] = toastMesh.Data();
    [surfpoints,edges,perm]=toastMesh.SurfData;
    vmcmesh.BH=[perm(edges(:,1)) perm(edges(:,2))];
end
