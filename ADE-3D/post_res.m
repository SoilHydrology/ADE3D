%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post processor - plot residuals
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

Rmax=max(abs(R_Omega),[],1);	% max |R_Omega|

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold;
  for k = Nt/10:Nt/10:Nt		% plot t=Nt*dt*1/10, Nt*dt*2/10, ...Nt*dt
  plot((1:Ne),R_Omega(:,k));
end	% k
xlabel('e');
ylabel('R_{\Omega}');
legend('Location','NorthEastOutside')
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
hold;
for e = 1:Ne
  plot((1:Nt),R_Omega(e,:));
end	% e
xlabel('k');
ylabel('R_{\Omega}');
legend('Location','NorthEastOutside')
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
plot(time, Rmax);
xlabel('time');
ylabel('Rmax = max(R_{\Omega})');
grid on;

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

Jmax=max(abs(J_Gamma),[],1);	% max |J_Gamma|

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
hold;
for k = Nt/10:Nt/10:Nt		% plot t=Nt*dt*1/10, Nt*dt*2/10, ...Nt*dt
  plot((1:Nf),J_Gamma(:,k));
end	% k
xlabel('f');
ylabel('J_{\Gamma}');
legend('Location','NorthEastOutside')
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
hold;
for f = 1:Nf
  plot((1:Nt),J_Gamma(f,:));
end	% f
xlabel('k');
ylabel('J_{\Gamma}');
legend('Location','NorthEastOutside')
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6)
plot(time, Jmax);
xlabel('time');
ylabel('Jmax = max(J_{\Gamma})');
grid on;

fclose('all');  % close all files

