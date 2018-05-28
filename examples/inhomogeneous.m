%% Creating an inhomogeneous medium: inhomogeneous.m
% This example demonstrates how to set up optical parameters that vary 
% as a function of position into a medium.  


%%

xsize =  10;	% width of the region [mm]
ysize =  10;	% height of the region [mm]
dh = 1;         % discretisation size [mm]
vmcmesh = createRectangularMesh(xsize, ysize, dh);

%% Give varying optical parameters
% The function 'createMedium' is a convenience function provided by 
% ValoMC that turns the fields in 'medium' into arrays. The sizes
% of the arrays are equal to the number of elements in the mesh, so that
% the optical parameters of each element can be set individually.
% 'findElementsInsideCircle' is an another convenience function provided
% to find the elements from the arrays.

% Set constant background coefficients
vmcmedium.absorption_coefficient = 0.01;     % absorption coefficient [1/mm]
vmcmedium.scattering_coefficient = 1.0;      % scattering coefficient [1/mm]
vmcmedium.scattering_anisotropy = 0.9;       % anisotropy parameter g of 
                                          % the Heneye-Greenstein scattering 
                                          % phase function [unitless]
vmcmedium.refractive_index = 1.3;            % refractive index [unitless]

% Create arrays of the coefficients
vmcmedium = createMedium(vmcmesh,vmcmedium);

% Find the elements that are inside a circle with a given radius and
% centre coordinates
radius = 2.5;                 % [mm]
centercoord = [0.0  0.0];     % [mm]
elements_of_the_circle = findElements(vmcmesh, 'circle', radius, centercoord);

% Assign a unique absorption coefficient to the circle
vmcmedium.absorption_coefficient(elements_of_the_circle) = 0.25;

%% A plot of the parameter distribution
% The resulting absorption coefficient distribution is shown in the figure
% below, along with the indices that are contained in "elements_of_the_circle".
% Shown is also a fictious circle with a radius of 2.5 mm.
% Note how some triangles intersect the circle, resulting in a poor
% representation of the circular shape. The representation can be improved
% by a better mesh or by incresing the discretisation size.
%
% <<circle.png>>

%% Create a light source
% Set a light source from boundary elements 4 to 7 with a cosinic directivity pattern 

vmcboundary.lightsource(4:7) = {'cosinic'};

%% Run the Monte Carlo simulation
% Use the parameters that were generated to run the simulation in the mesh

options.photon_count = 1e7;

solution = ValoMC(vmcmesh, vmcmedium, vmcboundary,options);

%% Plot the solution 
% Note how the area with the higher absorption coefficient
% affects the fluence distribution (cf. 'Simple Example')

hold on;
patch('Faces',vmcmesh.H,'Vertices',vmcmesh.r,'FaceVertexCData', solution.element_fluence, 'FaceColor', 'flat', 'LineWidth',1.5);
xlabel('[mm]');
ylabel('[mm]');
c = colorbar;                       % create a colorbar
c.Label.String = 'Fluence [J/mm^2]';
hold off

