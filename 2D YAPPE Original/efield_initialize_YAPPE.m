%This function resamples the input electric field (s.input.E_in) from the
%input radial and local time axes (s.input.r_in and s.input.xi_in) onto the
%radial and local time axes used by the code in the main computation loop
%(s.g.r and s.g.xi).

function [s] = efield_initialize_YAPPE(s)

%generate/load in linearly sampled input
if strcmp(s.input.infield.type, 'gauss')
    
    s = gauss_maker_YAPPE(s); %this makes a linearly sampled gaussian beam
    
    %interpolate the E_in onto the simulation xi and r grid points
    r_in_mat = bsxfun(@times,s.input.r_in,ones(size(s.input.xi_in)));
    xi_in_mat = bsxfun(@times,ones(size(s.input.r_in)),s.input.xi_in);
    s.f.E = interp2(xi_in_mat,r_in_mat,s.input.E_in,s.g.xi,s.g.r);
elseif strcmp(s.input.infield.type, 'custom')
    
    v = load(s.input.infield.path); %this loads in some time domain electric field
%     s.input.E_in = v.E_in;
%     s.input.r_in = v.r_in;
%     s.input.xi_in = v.xi_in;
    
    s.f.E = v.last_E_field;
    
end

% % delta_phi = 10*pi/s.input.r_pts;
% % for j = 1:s.input.r_pts
% %    s.f.E(j,:) = s.f.E(j,:)-1i*(j-1)*delta_phi;
% % end

%also initialize the spectral electric field
s.f.Ef = fft(s.f.E,[],2);
s.f.Ef = s.f.H*s.f.Ef;

end