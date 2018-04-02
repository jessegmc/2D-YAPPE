%This function creates a gaussian with a specified transverse size (waist),
%pulse length (tau), energy (energ) and quadratic focusing (f), and medium
%wavelength (kmed). If you don't want focusing, enter f= inf. The gaussian
%is defined on the grid r, xi,

function[s] = gauss_maker_YAPPE(s)

%create radial and temporal axes - (oversample to account for nonlinear spacing of radial axis)
s.input.xi_in = linspace(0, s.g.xi_range, 10*s.g.xi_pts); %local time axis in seconds
s.input.r_in = linspace(0, s.g.r_extent, 10*s.g.r_pts)'; %local radial axis
dxi = s.input.xi_in(2)-s.input.xi_in(1);
dr = s.input.r_in(2) - s.input.r_in(1);

%gaussian time constant
tau = s.input.infield.tfwhm /sqrt(2*log(2));

%create the desired gaussian shape
Er_env(:,1) = exp(-(s.input.r_in/s.input.infield.waist).^2).*exp(-1i*s.g.kcen*s.input.r_in.^2/(2*s.input.infield.f)); %changed s.g.k(1) to s.g.kcen, but f = inf so doesn't really matter
Et_env(1,:) = exp(-((s.input.xi_in - s.input.xi_in(end)/2)/tau).^2);
E_env = bsxfun(@times,Er_env,Et_env);

%normalize gaussian to create desired energy (assumes I = abs(E.^2) and integral of I = energy)
I = abs(E_env.^2);
nom_energ = 2*pi*dr*dxi*sum(sum(bsxfun(@times,I,s.input.r_in)));
s.input.E_in_env = sqrt(s.input.infield.energ/nom_energ)*E_env;

end