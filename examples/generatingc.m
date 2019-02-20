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

!../MC2D netgen_test_input.txt netgen_test_output.txt

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

