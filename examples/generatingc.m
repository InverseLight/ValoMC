%% Generating input for the external executable: generatingc.m
%
% Calculations with complex geometries are convenient to perform
% using the external executable that is compiled from C++ code
% (see installation instructions).  This example shows how the
% initial conditions can be set up using Matlab and stored to a text file
% for a parallelized calculation on a computer cluster, for example.  The code is
% identical to the Netgen example, but the computation is done using
% the external executable.
%
% *For this example to work the external executables must be compiled.*
%  *See the homepage/installation for instructions how to compile the external executales*

%%

clear all;

%% Perform simulation initialization as in the Netgen example
% Import the NetGen file
% Netgen meshes can be imported using 'importNetGenMesh'. In addition
% to the mesh structure, it returns the regions in the vol file as
% cell arrays. Each row in these arrays contains a vector that holds the
% indices of a region (in the medium or in the boundary). Netgen files
% may also contain names for the regions. These are returned in
% 'region_names' and 'boundary_names' and can be used to find the
% correct indices for each region (see 'Set optical parameters and light sources').

[vmcmesh regions region_names boundaries vmcboundary_names] = importNetGenMesh('square_with_two_circles.vol');

% Set optical parameters and light sources
% The return values can be used to assign optical coefficients, 
% lightsources and other conditions.

indices_for_background = cell2mat(regions(1));
indices_for_circles = cell2mat(regions(2));
indices_for_outer_boundary = [cell2mat(boundaries(2)); cell2mat(boundaries(1))];

vmcmesh.BH = vmcmesh.BH(indices_for_outer_boundary,:); 
indices_for_lightsource=1:size(cell2mat(boundaries(2)),1);

% In Matlab 2016b and later it is possible to find indices using
%
% indices_for_lightsource = cell2mat(boundaries(find(contains(boundary_names,'lightsource'))));
% indices_for_circles = cell2mat(regions(find(contains(region_names,'circles'))));
%
% i.e. strings can be used to extract regions.

vmcmedium.absorption_coefficient(indices_for_background) = 0.01;   % absorption coefficient [1/mm]
vmcmedium.scattering_coefficient(indices_for_background) = 1.3;    % scattering coefficient [1/mm]
vmcmedium.scattering_anisotropy(indices_for_background) = 0.9;     % scattering anisotropy parameter [unitless]
vmcmedium.refractive_index(indices_for_background) = 1.3;          % refractive index [unitless]

vmcmedium.absorption_coefficient(indices_for_circles) = 0.09;
vmcmedium.scattering_coefficient(indices_for_circles) = 1.3;
vmcmedium.scattering_anisotropy(indices_for_circles) = 0.5;
vmcmedium.refractive_index(indices_for_circles) = 1.5;

vmcboundary.lightsource(indices_for_lightsource) = {'cosinic'};    % cosine directed light profile

vmcoptions.photon_count = 2e6; % set the desired photon count

%% Save the Monte Carlo simulation input

% Export the input file 'netgen_test_input.txt'
exportValoMC('netgen_test_input.txt',vmcmesh, vmcmedium, vmcboundary, vmcoptions);

%% Run the external executable
% The input file is used to launch an external executable using the ! 
% operator in Matlab. Note that the calculation could be done on a computing 
% cluster aswell and no Matlab is needed.

% This assumes the c++ code has been compiled.
% In Windows, MC2D.a should be replaced with MC2D.exe

!./MC2D.a netgen_test_input.txt netgen_test_output.txt

%% Load the simulation output using importValoMC
% importValoMC can be used to retrieve the problem definition and 
% the simulation output from the external executable
[vmcmesh, vmcmedium, vmcboundary, options, solution] = importValoMC('netgen_test_input.txt', 'netgen_test_output.txt');

% Plot the solution

hold on;
patch('Faces',vmcmesh.H,'Vertices',vmcmesh.r,'FaceVertexCData', solution.element_fluence, 'FaceColor', 'flat', 'EdgeColor', 'none');
xlabel('[mm]');
ylabel('[mm]');
c = colorbar;                       
title('Fluence [W/mm^2]');
hold off

