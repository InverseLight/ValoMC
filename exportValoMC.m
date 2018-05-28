function exportValoMC(filename, vmcmesh, vmcmedium, vmcboundary, vmcoptions)
    vmcoptions.export_filename = filename;
    ValoMC(vmcmesh, vmcmedium, vmcboundary, vmcoptions);    
end
