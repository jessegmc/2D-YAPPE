% For use in a for loop. Define a global struct s constructed from the
% input deck, and this function will run YAPPE for the specified initial
% conditions. The struct s has to be global so that the matlab ODE function
% call can use information from s.

function[] = batch_runner_YAPPE()

%%%%%% INITIALIZATIONS %%%%%%%%%

global s
s = fund_const_YAPPE(s); %fundamental constants
s = medium_property_YAPPE(s); %medium parameters
s = grids_YAPPE(s); %make grids and interpolate input E-field onto new grids
s = precompute_YAPPE(s); %precompute useful arrays (Hankel matrix, Kz, Q, various coefficiencts...)
s = efield_initialize_YAPPE(s); %generate input E-fields (spatiotemporal and spectral)

%initialize outputs (helps with speed and acts as a check that there is enough memory for outputs)
s.f.E_env_out = zeros(s.input.r_pts,s.input.xi_pts,length(s.g.zout)); %initialize output E-field matrix
s.f.E_env_out(:,:,1) = s.f.E_env; %first entry is input field

%Keep track of Ef
s.f.Ef_out = zeros(s.input.r_pts,s.input.xi_pts,length(s.g.zout)); %initialize output Ef-field matrix
s.f.Ef_out(:,:,1) = s.f.Ef; 

s.f.rho_out_r = zeros(s.input.r_pts,length(s.g.zout)); %initialize output plasma matrix
s.f.rho_out_xi = zeros(s.input.xi_pts,length(s.g.zout)); %initialize output plasma matrix

s.count = 1; %used for monitoring # of ode45 function calls
opts = odeset('RelTol',s.input.RelTol,'AbsTol',s.input.AbsTol); %set solution tolerances

%%%%%% MAIN COMPUTATION LOOP %%%%%%%%%

for m = 2:length(s.g.zout)
    
    initial_Ef = s.f.Ef(:); %convert e-field matrix into vector for ode45 solve, at this point Ef = Af
    [~,x] = ode45( @(t,y) PNL_step_YAPPE(t,y), [0 s.input.outperiod], initial_Ef, opts); %solve ODE
    
    s.f.Af = reshape(x(end,:)', [s.input.r_pts, s.input.xi_pts]); %convert output of ODE solve to matrix, at this point Af=/=Ef
    clear x
    s.f.Af = conj(s.f.Af); %ODE45 has a "dagger" where they should have a "transpose", must conjugate
    s.f.Af_new = s.f.Af.*s.f.lin_prop; %apply linear propagator "realignment", now Af=Ef
    s.f.Ef = s.f.Af_new;
    
    %convert from spectral to spatiotemporal and write to output array
    s.f.Ef_out(:,:,m) = s.f.Ef;
    s.f.E_env_out_shift = ifft(s.f.Ef,[],2);
    s.f.E_env_out_shift = s.f.H*s.f.E_env_out_shift;
%     s.f.E_env_out(:,:,m) = s.f.E_env_out_shift.*exp(-1i*s.g.axis_shift*s.g.xi); %factor out unwanted frequency component
    s.f.E_env_out(:,:,m) = s.f.E_env_out_shift.*exp(-1i*s.g.num_pts_to_shift/s.g.dxi*2*pi/s.g.xi_pts*s.g.xi); %factor out unwanted frequency component
    %write out the the local time-integrated plasma
    s.f.rho_out_r(:,m) = squeeze(sum(s.f.rho*s.g.dxi,2)); %radially resolved time integrated plasma density
    s.f.rho_out_xi(:,m) = squeeze(s.f.rho(1,:)); %axial lineout of plasma density
    
    disp(strcat( 'propagated to z =', 32, num2str(s.g.zout(m)), 32, 'cm'))
    
end

%%%%%% WRITING OUT %%%%%%%%%

%save the outputs configuration
outputnam = strcat(s.input.outpath,'full_output.mat');
if exist(s.input.outpath,'file')==0
    mkdir(s.input.outpath);
end
save(outputnam, 's');

disp('a million coupled ODEs cried out... and were solved')
end
