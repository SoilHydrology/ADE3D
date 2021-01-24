%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post processor - arrange FEM-3D in matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all 
close all;	% close all the open figure windows

%%%%%%%%%%% read meshinp.txt %%%%%%%%%%%%
fid = fopen('meshinp.txt','r');
line = fgets(fid);	% read line (Darcy)
line = fgets(fid);	% read line (Xmin, Xmax, Ymin, Ymax, Zmin, Zmax)
line = fgets(fid);	% read line (Nx, Ny, Nz, Ns)
v = sscanf(line, '%i');
Nx = v(1);
Ny = v(2);
Nz = v(3);

%%%%%%%%%%% read finp.txt %%%%%%%%%%%%
fid = fopen('finp.txt','r');

line = fgets(fid);	% read line (Ne, Nn, Nb, Nm, Np)
line = fgets(fid);	% read line (Nd, Ng, Ns)
v = sscanf(line, '%i');
Ns = v(3);
line = fgets(fid);	% read line (Kappa, tc, por)
line = fgets(fid);	% read line (tmax, dt, dto)
v = sscanf(line, '%e');
tmax = v(1);
dto = v(3);
Nto = (tmax + 1e-6) / dto;

%%%%%%%%%% read post.msh %%%%%%%%%%%
fid = fopen('post.msh','r');
for m =1:6
  line = fgets(fid);	% read 6 lines
end	% m

% read nodal coordinates
for k = 1:Nz+1
  for j = 1:Ny+1
    for i = 1:Nx+1
      line = fgets(fid);	% read line (n, x, y, z)
%      v = sscanf(line, '%i %e %e %e');
    end	% i
  end	% j
end	% k

for m =1:5
  line = fgets(fid);	% read 5 lines
end	% m

% read elements
for k = 1:Nz
  for j = 1:Ny
    for i = 1:Nx
      line = fgets(fid);	% read line (elm-number elm-type number-of-tags < tag > ... node-number-list)
     end	% i
  end	% j
end	% k

for m =1:12
  line = fgets(fid);	% read 12 lines
end	% m

% read elements data (V)
for k = 1:Nz
  for j = 1:Ny
    for i = 1:Nx
      line = fgets(fid);	% read line (elm-number number-of-nodes-per-element value ... Vx, Vy, Vz)
    end	% i
  end	% j
end	% j

for m =1:12
  line = fgets(fid);	% read 12 lines
end	% m

% read elements data (D)
for k = 1:Nz
  for j = 1:Ny
    for i = 1:Nx
      line = fgets(fid);	% read line (elm-number number-of-nodes-per-element value ... Dxx, Dxy,.. Dzz)
    end	% i
  end	% j
end	% k

for m =1:12
  line = fgets(fid);	% read 12 lines
end	% m

% time loop
for t=0:Nto

% read nodal data (Cs)
  for s = 1:Ns
    for k = 1:Nz+1
      for j = 1:Ny+1
        for i = 1:Nx+1
          line = fgets(fid);	% read line (node-number value ...)
          v = sscanf(line, '%i %e');
          C(t+1,s,k,j,i) = v(2);
        end	% i
      end	% j
    end		% k

    for m =1:12
      line = fgets(fid);	% read 12 lines
    end		% m
  end		% s

% read nodal data (Ts)
  for s = 1:Ns
    for k = 1:Nz+1
      for j = 1:Ny+1
        for i = 1:Nx+1
          line = fgets(fid);	% read line (node-number value ...)
          v = sscanf(line, '%i %e');
          T(t+1,s,k,j,i) = v(2);
        end	% i
      end	% j
    end		% k

    for m =1:12
      line = fgets(fid);	% read 12 lines
    end		% m
  end		% s
end		% t

% % read nodal data (q_i)
%  for k = 1:Nz+1
%   for j = 1:Ny+1
%     for i = 1:Nx+1
%       line = fgets(fid);	% read line (node-number value ...)
%       v = sscanf(line, '%i %e %e %e');
%       qx(j,i)= v(2);
%       qy(j,i)= v(3);
%       qz(j,i)= v(4);
%     end	% i
%   end		% j
% end		% k

% C1 = reshape(C(:,1,:,:,:),[(Nz+1)*(Ny+1)*(Nx+1),t+1,]);
% C2 = reshape(C(:,2,:,:,:),[(Nz+1)*(Ny+1)*(Nx+1),t+1,]);
% C3 = reshape(C(:,3,:,:,:),[(Nz+1)*(Ny+1)*(Nx+1),t+1,]);

fclose('all');  % close all files
