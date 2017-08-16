%loads material parameters into substruct s.mat.

function [s] = medium_property_YAPPE(s)

if strcmp(s.input.medium, 'air')
    
    s.mat.tau_c = 3.5e-13; %electron collision time in seconds
    
    s.mat.N0 = 2.504e19; %number density of air in cm^-3
    s.mat.NO_N2 = 2e19; %number density of N2 in the air
    s.mat.NO_O2 = 5e18; %number density of O2 in the air
    
    s.mat.Eg_N2 = 2.496e-18; %ionization energy in J of N2 (15.58 eV)
    s.mat.Eg_O2 = 1.934e-18; %ionization energy in J of O2 (12.07 eV)
    
    s.mat.R_N2 = 2.5e4; %Rate prefactor (experimentally fitted)
    s.mat.R_O2 = 2.8e6; %Rate prefactor (experimentally fitted)
    
    s.mat.IT = 1e13; %Intensity factor (experimentally fitted)
    
    %     s.mat.recomb = 5e-13*1e6; %recombination rate in cm^3/s (was "s.a3")
    
    s.mat.alpha_N2 = 7.5; %Number of photons to ionize N2
    s.mat.alpha_O2 = 6.5; %Number of photons to ionize O2
    
    s.mat.n2_N2 = 7.9e-20; %nonlinear index of N2 - valid at 800nm in cm^2/W
    s.mat.n2_O2 = 8.1e-20; %nonlinear index of O2 - valid at 800nm in cm^2/W
    s.mat.n2 = 0.78*s.mat.n2_N2+0.21*s.mat.n2_O2; %rough estimate of average n2 %7.8e-20
    
end

end