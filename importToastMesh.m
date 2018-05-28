function vmcmesh = importToastMesh(toastMesh)
   
    if not(isa(toastMesh, 'toastMesh'))
        error('Could not recognize toast mesh');
    end
    element=toastMesh.Element(1);
    if(element.elid ~= 1)
        error('Element type currently not supported');
    end
    [vmcmesh.r, vmcmesh.H] = toastMesh.Data();
    [surfpoints,edges,perm]=toastMesh.SurfData;
    vmcmesh.BH=[perm(edges(:,1)) perm(edges(:,2))];
end 
