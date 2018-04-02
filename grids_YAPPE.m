%create grids and grid parameters and load into substruct s.g

function [s] = grids_YAPPE(s)

%lab frame axis properties
s.g.z_range = s.input.z_extent; %propagation range in cm
s.g.zout = 0:s.input.outperiod:s.g.z_range; %output propagation locations in cm

if ~s.input.userSuppliedFreqAxis %user supplies time extent
    %local time axis grid and grid properties
    s.g.xi_range = s.input.xi_extent; %temporal range in seconds
    s.g.xi_pts = s.input.xi_pts; %number of grid pts along local time axis
    s.g.xi = linspace(0, s.g.xi_range, s.input.xi_pts); %local time axis in seconds
    s.g.dxi = s.g.xi(2) - s.g.xi(1); %local time axis grid spacing in seconds
    
    %spectral axis grid and grid properties
    s.g.omg_cen = 2*pi*s.cgs.c/s.input.lambda_vac; %central angular frequency in rad/s
    s.g.domg = 2*pi/(s.g.xi_range); %angular frequency spacing in rad/s
    s.g.omg = s.g.domg*((1:s.g.xi_pts) - round((s.g.xi_pts+1)/2)); %angular frequency axis in rad/s
    s.g.omg = s.g.omg + s.g.omg_cen; %center axis about central frequency
    s.g.omg = ifftshift(s.g.omg,2); %shift the axis
    s.g.omg_cen_loc = 1; %location of central frequency of the pulse within the frequency axis
    s.g.axis_shift = 0; %difference between central frequencies of the pulse and of the axis
    s.g.num_pts_to_shift = 0; %number of points to shift the frequency axis by
    
else %user supplies angular frequency range
    s.g.xi_pts = s.input.xi_pts; %number of grid pts along local time axis
    s.g.omg_cen = 2*pi*s.cgs.c/s.input.lambda_vac; %central angular frequency of pulse in rad/s
    
    s.g.bandwidth = s.input.omg_max-s.input.omg_min; %angular frequency bandwidth
    s.g.domg = s.g.bandwidth/(s.g.xi_pts-1); %angular frequency spacing in rad/s -- the minus 1 accounts for min and max inclusion in the next line
    s.g.omg = s.input.omg_min:s.g.domg:s.input.omg_max; %angular frequency axis in rad/s
    s.g.omg_cen_loc = round((s.g.omg_cen-s.input.omg_min)/s.g.domg)+1; %the location of omg_cen within the omg axis
    
    s.g.omg_cen_test = s.g.omg(s.g.omg_cen_loc); %test to see if omg_cen is a point on omg axis
    s.g.omg = s.g.omg + (s.g.omg_cen-s.g.omg_cen_test); %shift entire omg axis so that omg_cen is always a point on the axis

    s.g.omg_cen_axis = (s.input.omg_max+s.input.omg_min+s.g.domg)/2; %central angular frequency of freq axis, must be even number of xi_pts
    s.g.axis_shift = s.g.omg_cen-s.g.omg_cen_axis; %angular frequency difference between the central frequency of the pulse and central frequency of the frequency axis
    
    s.g.num_pts_to_shift = round(s.g.axis_shift/s.g.domg); %Number of places to shift freq axis 
    
    s.g.omg = ifftshift(s.g.omg,2); %shift the axis
    if(s.g.omg_cen_loc<=round(s.g.xi_pts/2))
        s.g.omg_cen_loc = s.g.omg_cen_loc+round(s.g.xi_pts/2); %account for change in location of omg_cen due to ifftshift
    else
        s.g.omg_cen_loc = s.g.omg_cen_loc-round(s.g.xi_pts/2); %account for change in location of omg_cen due to ifftshift
    end
    
    %local time axis grid and grid properties
    s.g.xi_range = 2*pi/s.g.domg; %temporal range in seconds
    s.g.xi = linspace(0, s.g.xi_range, s.input.xi_pts); %local time axis in seconds
    s.g.dxi = s.g.xi(2) - s.g.xi(1); %local time axis grid spacing in seconds
end

%wavenumbers
s.g.kvac = s.g.omg/s.cgs.c; %vaccuum wavenumbers
s = dispersion_YAPPE(s); %get permittivity
s.g.n = sqrt(s.g.perm); %linear index
s.g.n0 = s.g.n(s.g.omg_cen_loc); %central index
s.g.k = s.g.omg.*s.g.n/s.cgs.c; %medium wavenumbers
s.g.kcen = s.g.k(s.g.omg_cen_loc); %central medium wavenumber

%group velocity
if(s.g.omg_cen_loc==1) %accounts for the case in which the central ang freq is the first frequency in the omega axis
    s.g.vg = ((s.g.k(2) - s.g.k(end))/(2*s.g.domg) )^(-1);
elseif(s.g.omg_cen_loc==s.g.xi_pts) %accounts for the case in which the central ang freq is the last frequency in the omega axis
    s.g.vg = ((s.g.k(1) - s.g.k(end-1))/(2*s.g.domg) )^(-1);
else %general case
    s.g.vg = ((s.g.k(s.g.omg_cen_loc+1) - s.g.k(s.g.omg_cen_loc-1))/(2*s.g.domg) )^(-1);
end

%transverse spatial grid and grid properties
s.g.r_pts = s.input.r_pts; %radial window length in cm
s.g.r_extent = s.input.r_extent;
s = hankel_sample_YAPPE(s); %this generates the bessel zero-spaced radial grid

end