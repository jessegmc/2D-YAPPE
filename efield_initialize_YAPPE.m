%This function resamples the input electric field (s.input.E_in) from the
%input radial and local time axes (s.input.r_in and s.input.xi_in) onto the
%radial and local time axes used by the code in the main computation loop
%(s.g.r and s.g.xi).

function [s] = efield_initialize_YAPPE(s)

%generate/load in linearly sampled input
if s.input.infield.type == 'gauss'
    
    s = gauss_maker_YAPPE(s); %this makes a linearly sampled gaussian beam
    
elseif s.input.infield.type == 'custom'
    
    v = load(s.input.infield.path); %this loads in some time domain electric field
    s.input.E_in = v.E_in;
    s.input.r_in = v.r_in;
    s.input.xi_in = v.xi_in;
    
end

%interpolate the E_in onto the simulation xi and r grid points
r_in_mat = bsxfun(@times,s.input.r_in,ones(size(s.input.xi_in)));
xi_in_mat = bsxfun(@times,ones(size(s.input.r_in)),s.input.xi_in);
s.f.E_env = interp2(xi_in_mat,r_in_mat,s.input.E_in_env,s.g.xi,s.g.r);

% s.f.E_env_freq_shift = s.f.E_env.*exp(1i*s.g.axis_shift*s.g.xi);
s.f.E_env_freq_shift = s.f.E_env.*exp(1i*s.g.num_pts_to_shift/s.g.dxi*2*pi/s.g.xi_pts*s.g.xi);

%If the central frequency of the frequency axis is not the same as the 
%central frequency of the pulse, we must multiply the time domain envelope 
%by the difference in frequenceies between pulse central frequency and 
%axis central frequency. This ensures that when we take the DFT, the new
%variable is centered around the central frequency of the pulse.

%also initialize the spectral electric field
s.f.Ef = fft(s.f.E_env_freq_shift,[],2);
s.f.Ef = s.f.H*s.f.Ef; %The DFT of the envelope of E gives us Ef when viewed through our frequency axis

% outputnam = strcat(s.input.outpath,'full_output.mat');
% if exist(s.input.outpath,'file')==0
%     mkdir(s.input.outpath);
% end
% save(outputnam, 's');

end