%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post processor - residuals RMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all 
close all

%%%%%%%%%%% read finp.txt %%%%%%%%%%%%
fid = fopen('finp.txt','r');

line = fgets(fid);	% read line (Ne, Nn, Nb, Nm, Np)
v = sscanf(line, '%i');
Ne = v(1);
Nn = v(2);
line = fgets(fid);	% read line (Nd, Ns, Npol, Nf)
v = sscanf(line, '%i');
Npol = v(3);
Nf   = v(4);
line = fgets(fid);	% read line (Sx, Kappa, tc, por)
line = fgets(fid);	% read line (tmax, dt, dto)
v = sscanf(line, '%e');
tmax = v(1);
dt   = v(2);
dto  = v(3);
Nt = round(tmax/dt );    % round to nearest integer
Nto= round(tmax/dto);    % round to nearest integer
line = fgets(fid);	% read line (theta)
line = fgets(fid);	% read line (ipar(2:3), ipar(6))
line = fgets(fid);	% read line (fpar(1:2), fpar(11))

%%%%%%%%%%% read resi.txt %%%%%%%%%%%%
R_Omega = zeros(Ne, Nt);
time   = zeros( 1, Nt);
fid = fopen('resi.txt','r');
for k = 1:Nt
  d = fscanf(fid, '%s %s',2);	% dummy
  time(k) = fscanf(fid, '%e',1);	% time
  d = fscanf(fid, '%s',1);	% dummy (header)
  for e = 1:Ne
    R_Omega(e, k) = fscanf(fid, '%e' ,1);	% R_Omega
  end	% e
end	% k

R_Omega_RMS = rms(R_Omega,2);	% RMS R_Omega for each element

figure(1);plot((1:Ne), R_Omega_RMS);xlabel('e');ylabel('R_{\Omega}^{RMS}');grid on;

%%%%%%%%%%% read jump.txt %%%%%%%%%%%%
J_Gamma = zeros(Nf, Nt);
fid = fopen('jump.txt','r');
for k = 1:Nt
  d = fscanf(fid, '%s %s',2);	% dummy
  time(k) = fscanf(fid, '%e',1);	% time
  d = fscanf(fid, '%s',1);	% dummy (header)
  for f = 1:Nf
    J_Gamma(f, k) = fscanf(fid, '%e' ,1);	% J_Gamma
  end	% f
end	% k

J_Gamma_RMS = rms(J_Gamma,2);	% RMS J_Gamma for each face

figure(2);plot((1:Nf), J_Gamma_RMS);xlabel('f');ylabel('J_{\Gamma}^{RMS}');grid on;

%%%%%%%%%%% read porder.txt %%%%%%%%%%%%
fid = fopen('porder.txt','r');
epol = zeros(1,Ne);
for e = 1:Ne
  line = fgets(fid);	% read line (e, epol(e))
  v = sscanf(line, '%i');
  epol(e) = v(2);
end	% e

figure(3);plot((1:Ne), epol);xlabel('e');ylabel('epol');grid on;

fclose('all');  % close all files

