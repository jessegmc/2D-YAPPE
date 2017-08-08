%This function creates the input struct for YAPPE.

function [s] = input_deck_YAPPE()

%grid parameters
s.input.xi_extent = 600*1e-15; %axial window length in seconds
s.input.xi_pts = 300; %number of pts along axial direction
% % s.input.r_extent = .2; %radial window length in cm
s.input.r_extent = .2*1e-2; %radial window length in m
s.input.r_pts = 600; %number of pts along radial direction
% % s.input.z_extent = 0.05; %extent of propagation in cm
s.input.z_extent = 0.5*1e-2; %extent of propagation in m

%load E-field options
% % s.input.lambda_vac = 800e-7; %vacuum central wavelength in cm
s.input.lambda_vac = 800e-9; %vacuum central wavelength in m
s.input.infield.type = 'gauss'; %choose between 'gauss' and 'custom'
s.input.infield.path = '~/matlab_scripts/hankel_FILA/E_in_example.mat'; %location of custom field
% % s.input.infield.waist = 40e-4; %gaussian waist in cm (used in gauss)
s.input.infield.waist = 40e-4*1e-2; %gaussian waist in m (used in gauss)
s.input.infield.tfwhm = 40e-15; %intensity fwhm in s (used in gauss)
s.input.infield.energ = .5e-6*(0.5*3e8*1.3282*8.854e-12); %beam energy in J (used in gauss)
% % s.input.infield.f = inf; %focusing length in cm (used in gauss)
s.input.infield.f = inf; %focusing length in m (used in gauss)

%choose propagation medium
s.input.medium = 'water';

%specify output path and output period
s.input.outpath = 'R:\Jesse GPU\MATLAB\YAPPE Output\2D\2D YAPPE mks\Run4\'; %output path
% % s.input.outperiod = .05; %output period in cm
s.input.outperiod = .05*1e-2; %output period in m

%toggle dispersion, plasma and n2 propagation modules
s.input.dispersion = 1;
s.input.plasma = 0;
s.input.n2 = 0;

%solution tolerances for ODE calls
s.input.RelTol = 1e-3;
s.input.AbsTol = 1e-6;

%absorbing boundaries
s.input.freqbd = 1; %toggle the absorbing frequency boundary on and off
s.input.freqbd_length = .005; %set absorption length in cm
s.input.freqbd_width = .1; %fractional value of boundary width (boundary width/total bandwidth)


end