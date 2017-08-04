%This function steps the solution forward w.r.t. the nonlinear
%polarizability. This is where almost all the computation time is. The load
%is fairly equally spread between 2 fft calls, 2 hankel transforms, the
%plasma integration (vector sums), and the 2 element-wise matrix
%exponentiation calls. Declaring the struct "s" as global does not have a
%performance cost. I observe roughly equal call times between the above
%operations for a 300x300 grid.

function[dEf_dz_vec] = PNL_step_YAPPE(z,Ef_vec)

global s

%reshape Ef back into a matrix
Ef = reshape(Ef_vec,s.input.r_pts,s.input.xi_pts);
Ef = Ef.*exp(1i*s.f.Kz_move*z);

%convert to spatiotemporal domain
E = ifft(Ef,[],2);
E = s.f.H*E;

%calculate intensity envelope
s.f.I = abs(E).^2;

%calculate plasma density in #/cm^3
s.f.rho = zeros(size(s.f.I));

%Plasma density is calculated from MPI with a prefactor and scaling term
%which accounts for tunneling and MPI in the 10-100 TW/cm^2 regime
%The model is fitted so N2 absorbs 7.5 photons and O2 absorbs 6.5 photons.

if s.input.plasma == 1 %this toggles the plasma module
    for m = 2:size(s.f.rho,2)
        rhoN2 = s.mat.NO_N2*s.mat.R_N2*(s.f.I(:,m-1)/s.mat.IT).^s.mat.alpha_N2;
        rhoO2 = s.mat.NO_O2*s.mat.R_O2*(s.f.I(:,m-1)/s.mat.IT).^s.mat.alpha_O2;
        s.f.rho(:,m) = s.f.rho(:,m-1) +s.g.dxi*(rhoN2+rhoO2);
    end
end

%calculate nonlinear susceptibility
%The NL susceptibility does not have a term that accounts for pulse energy 
%loss due to plasma formation
s.f.chiNL = s.input.n2*( s.NL.b1*s.f.I ) + s.input.plasma*( s.NL.b2*s.f.rho );

%calculate nonlinear polarizability
s.f.PNL = s.f.chiNL.*E;

%convert nonlinear polarizability to spectral domain
s.f.PNLf = fft(s.f.PNL,[],2);
s.f.PNLf = s.f.H*s.f.PNLf;

%z derivative in matrix form
dEf_dz =  0.5*1i*s.f.Q.*s.f.PNLf.*exp(-1i*s.f.Kz_move*z);

%z derivative in vector form
dEf_dz_vec = dEf_dz(:);
s.count = s.count+1;

end
