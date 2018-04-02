%This function steps the solution forward w.r.t. the nonlinear
%polarizability. This is where almost all the computation time is. The load
%is fairly equally spread between 2 fft calls, 2 hankel transforms, the
%plasma integration (vector sums), and the 2 element-wise matrix
%exponentiation calls. Declaring the struct "s" as global does not have a
%performance cost. I observe roughly equal call times between the above
%operations for a 300x300 grid.

function[dAf_dz_vec] = PNL_step_YAPPE(z,Ef_vec)

global s

%reshape Ef back into a matrix
Ef = reshape(Ef_vec,s.input.r_pts,s.input.xi_pts);
Ef = Ef.*exp(1i*s.f.Kz_move*z); %Apply linear propagator to shift z position of the pulse

%convert to spatiotemporal domain
E_env_freq_shift = ifft(Ef,[],2);
E_env_freq_shift = s.f.H*E_env_freq_shift;

% E_env = E_env_freq_shift.*exp(-1i*s.g.axis_shift*s.g.xi); %taking the DFT of Ef returns E_env multiplied by exp(1i*s.g.axis_shift*s.g.xi), must factor out the exponential to single out E_env 
E_env = E_env_freq_shift.*exp(-1i*s.g.num_pts_to_shift/s.g.dxi*2*pi/s.g.xi_pts*s.g.xi);  %taking the DFT of Ef returns E_env multiplied by exp(1i*s.g.axis_shift*s.g.xi), must factor out the exponential to single out E_env 


% %save the outputs configuration
% outputnam1 = strcat(s.input.outpath,'E_env.mat');
% outputnam2 = strcat(s.input.outpath,'Ef.mat');
% outputnam3 = strcat(s.input.outpath,'full_output.mat');
% 
% if exist(s.input.outpath,'file')==0
%     mkdir(s.input.outpath);
% end
% save(outputnam1, 'E_env');
% save(outputnam2, 'Ef');
% save(outputnam3, 's');


%calculate intensity envelope
s.f.I = abs(E_env).^2;
Ipow = s.f.I.^(s.mat.pow);

%calculate plasma density in #/cm^3
s.f.rho = zeros(size(s.f.I));

if s.input.plasma == 1 %this toggles the plasma module
    for m = 2:size(s.f.rho,2)
        s.f.rho(:,m) = s.f.rho(:,m-1) + s.g.dxi*( s.ion.a1*s.f.I(:,m-1).*s.f.rho(:,m-1) + s.ion.a2*Ipow(:,m-1) + s.ion.a3*s.f.rho(:,m-1).^2 );
    end
end

%calculate nonlinear susceptibility
s.f.chiNL = s.input.n2*( s.NL.b1*s.f.I ) + s.input.plasma*( s.NL.b2*s.f.rho + s.NL.b3*s.f.I.^(s.mat.pow-1) );

%calculate nonlinear polarizability
s.f.PNL_env = s.f.chiNL.*E_env;

% s.f.PNL_freq_shift = s.f.PNL_env.*exp(1i*s.g.axis_shift*s.g.xi); 
s.f.PNL_freq_shift = s.f.PNL_env.*exp(1i*s.g.num_pts_to_shift/s.g.dxi*2*pi/s.g.xi_pts*s.g.xi); 
%must multiply spatiotemporal domain polarizability by exp(1i*s.g.axis_shift*s.g.xi) so when we 
%take the DFT, the frequency domain is shifted to align with our frequency axis

%convert nonlinear polarizability to spectral domain
s.f.PNLf = fft(s.f.PNL_freq_shift,[],2);
s.f.PNLf = s.f.H*s.f.PNLf; 

%z derivative in matrix form
dAf_dz =  0.5*1i*s.f.Q.*s.f.PNLf.*exp(-1i*s.f.Kz_move*z);

%z derivative in vector form
dAf_dz_vec = dAf_dz(:);
s.count = s.count+1;

end
