%% Working with NetGen: netgentest.m
%
% This example demonstrates how to import a mesh from Netgen. The python
% source code for Netgen that generates the mesh can found in the
% examples/netgen_square_with_two_circles.py. The Python source code can be
% viewed <netgen_square_with_two_circles.py here>

%% Import the NetGen file
% Netgen meshes can be imported using 'importNetGenMesh'. In addition to
% the mesh structure, it returns the regions and boundaries in the vol file
% as cell arrays. If the second argument is set to 'false', a new boundary
% will be generated and the one in the file will not be used. It is
% recommended since the original boundary oftain contains boundary elements
% that are between normal elements. This is currently not supported in
% ValoMC.

clear all;

[vmcmesh regions region_names boundaries boundary_names] = importNetGenMesh('square_with_two_circles.vol', false);


%% Find indices

indices_for_background = cell2mat(regions(find(strcmp(region_names,'background'))));
indices_for_circles = cell2mat(regions(find(strcmp(region_names,'circles'))));
indices_for_lightsource = cell2mat(boundaries(find(strcmp(boundary_names,'lightsource'))));


%% Set optical parameters and light sources using the indices

vmcmedium.absorption_coefficient(indices_for_background) = 0.01;   % absorption coefficient [1/mm]
vmcmedium.scattering_coefficient(indices_for_background) = 1.3;    % scattering coefficient [1/mm]
vmcmedium.scattering_anisotropy(indices_for_background) = 0.9;     % scattering anisotropy parameter [unitless]
vmcmedium.refractive_index(indices_for_background) = 1.3;          % refractive index [unitless]

vmcmedium.absorption_coefficient(indices_for_circles) = 0.09;
vmcmedium.scattering_coefficient(indices_for_circles) = 1.3;
vmcmedium.scattering_anisotropy(indices_for_circles) = 0.5;
vmcmedium.refractive_index(indices_for_circles) = 1.5;

vmcboundary.lightsource(indices_for_lightsource) = {'cosinic'};   


%% Run the Monte Carlo simulation
solution = ValoMC(vmcmesh, vmcmedium, vmcboundary);

%% Plot the solution

figure('rend','painters','pos',[10 10 1200 400])

h = subplot(1,2,1);
hold on;
patch('Faces',vmcmesh.H,'Vertices',vmcmesh.r,'FaceVertexCData', vmcmedium.absorption_coefficient(:), 'FaceColor', 'flat', 'EdgeColor','none');
xlabel('[mm]');
ylabel('[mm]');
c = colorbar;                       
hold off
title('Absorption coefficient [1/mm]');

h=subplot(1,2,2);
hold on;
patch('Faces',vmcmesh.H,'Vertices',vmcmesh.r,'FaceVertexCData', solution.element_fluence, 'FaceColor', 'flat', 'EdgeColor', 'none');
xlabel('[mm]');
ylabel('[mm]');
c = colorbar;                       
title('Fluence [W/mm^2]');
hold off

