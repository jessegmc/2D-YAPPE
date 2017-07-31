%This function steps the solution forward w.r.t. the nonlinear
%polarizability. This is where almost all the computation time is. The load
%is fairly equally spread between 2 fft calls, 2 hankel transforms, the
%plasma integration (vector sums), and the 2 element-wise matrix
%exponentiation calls. Declaring the struct "s" as global does not have a
%performance cost. I observe roughly equal call times between the above
%operations for a 300x300 grid.

function[dEf_dz_vec] = PNL_step_YAPPE_single(z,Ef)

global s

%reshape Ef back into a matrix
% %   Ef = reshape(Ef_vec,s.input.r_pts,s.input.xi_pts);
Ef = reshape(Ef,s.input.x_pts,s.input.y_pts,s.input.xi_pts);
% s.f.Kz_move = single(s.f.Kz_move(:));
Ef = single(Ef);
Ef = Ef.*exp(1i*s.f.Kz_move*z);

%convert to spatiotemporal domain
% %     E = ifft(Ef,[],2);
% %     E = s.f.H*E;

Ef = ifftn(Ef);
% a = sum(sum(sum(isnan(Ef))))
%calculate intensity envelope
s.f.I = abs(Ef).^2;
% b = sum(sum(sum(isnan(s.f.I))))
% b2 = sum(sum(sum(isinf(s.f.I))))

Ipow = double(s.f.I).^(s.mat.pow);
% c = sum(sum(sum(isnan(Ipow))))
% c2 = sum(sum(sum(isinf(Ipow))))
s.f.I = gather(s.f.I);
Ipow = gather(Ipow);
%calculate plasma density in #/cm^3
s.f.rho = zeros(size(s.f.I),'single');

if s.input.plasma == 1 %this toggles the plasma module
    % %         for m = 2:size(s.f.rho,2)
    % %             s.f.rho(:,m) = s.f.rho(:,m-1) + s.g.dxi*( s.ion.a1*s.f.I(:,m-1).*s.f.rho(:,m-1) + s.ion.a2*Ipow(:,m-1) + s.ion.a3*s.f.rho(:,m-1).^2 );
    % %         end
    for m = 2:size(s.f.rho,3)
        s.f.rho(:,:,m) = s.f.rho(:,:,m-1) + s.g.dxi*( s.ion.a1*s.f.I(:,:,m-1).*s.f.rho(:,:,m-1) + s.ion.a2*Ipow(:,:,m-1) + s.ion.a3*s.f.rho(:,:,m-1).^2 );
    end
end
clear Ipow;
s.f.I = gpuArray(s.f.I);
% d = sum(sum(sum(isnan(s.f.rho))))
% d2 = sum(sum(sum(isinf(s.f.rho))))
%calculate nonlinear susceptibility
s.f.chiNL = s.input.n2*( s.NL.b1*double(s.f.I) ) + s.input.plasma*( s.NL.b2*double(s.f.rho) + s.NL.b3*double(s.f.I).^(s.mat.pow-1) );

% e = sum(sum(sum(isnan(s.f.chiNL))))
% e2 = sum(sum(sum(isinf(s.f.chiNL))))
s.f = rmfield(s.f,'rho');

%calculate nonlinear polarizability
s.f.PNL = s.f.chiNL.*Ef;
% f = sum(sum(sum(isnan(s.f.PNL))))
% f2 = sum(sum(sum(isinf(s.f.PNL))))
s.f = rmfield(s.f,'chiNL');

%convert nonlinear polarizability to spectral domain
% %     s.f.PNLf = fft(s.f.PNL,[],2);
% %     s.f.PNLf = s.f.H*s.f.PNLf;

s.f.PNL = fftn(s.f.PNL);
% g = sum(sum(sum(isnan(s.f.PNL))))
%     clear s.f.PNL;

%z derivative in matrix form
dEf_dz_vec =  0.5*1i*s.f.Q.*s.f.PNL.*exp(-1i*s.f.Kz_move*z);
% h = sum(sum(sum(isnan(dEf_dz_vec))))
s.f = rmfield(s.f,'PNL');
%z derivative in vector form
dEf_dz_vec = dEf_dz_vec(:);
s.count = s.count+1;

end
