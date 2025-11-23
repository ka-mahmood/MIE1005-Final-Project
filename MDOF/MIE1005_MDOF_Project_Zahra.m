 clear; clc;

%%  1. PEMFC INTERNAL 7×7 MODEL (UNCOUPLED)

M = diag(0.001*[806.4,(22.275-0.495),(2.937+0.4185),0.165,...
               (0.145+0.4185),(22.275-0.495),806.4]);

% Material stiffness
k_endplate   = 4.788e10; 
k_BPP        = 5.775e10;
k_rib        = 1.115e10;
k_edge       = 5.775e10;
k_GDLa       = 9.6525e11;
k_GDLc       = 8.52e10;
k_PEM        = 7.953e9;
k_clamp      = 6.015e9;
k_gasket     = 3.056e10;
n = 1;

k_mount_top = 1e6;
k_endplate = (k_endplate * k_mount_top)/(k_endplate + k_mount_top);

k_plate     = ((1/(k_edge*2+k_rib))+(1/k_BPP))^-1;           
k_bpp_gdlc  = 2*k_gasket + n*k_GDLc;          
k_bpp_gdla  = 2*k_gasket + n*k_GDLa;       
k_GDL_pem   = (2+n)*k_PEM;

% stiffness matrix
K = [k_plate+k_clamp  -k_plate           0           0           0           0     -k_clamp;
    -k_plate          k_plate+k_bpp_gdlc -k_bpp_gdlc 0           0           0      0;
     0               -k_bpp_gdlc        k_bpp_gdlc+k_GDL_pem  -k_GDL_pem     0      0      0;
     0                0                -k_GDL_pem   k_GDL_pem+k_GDL_pem -k_GDL_pem  0      0;
     0                0                 0          -k_GDL_pem  k_GDL_pem+k_bpp_gdla -k_bpp_gdla 0;
     0                0                 0           0          -k_bpp_gdla  k_bpp_gdla+k_plate -k_plate;
    -k_clamp          0                 0           0           0          -k_plate  k_plate+k_clamp];

%% ADD VEHICLE SUSPENSION → EXTENDED 9×9 MODEL

m_s  = 2500;
m_us = 200;

k_s  = 80000;
c_s  = 350;
k_t  = 0.5e6;
c_t  = 15020;

M9 = blkdiag(M, m_s, m_us);
K9 = blkdiag(K, zeros(2));

% mount connection
K9(1,1) = K9(1,1) + k_endplate;
K9(7,7) = K9(7,7) + k_endplate;
K9(7,8) = K9(7,8) - k_endplate;   K9(8,7) = K9(8,7) - k_endplate;
K9(8,8) = K9(8,8) + k_endplate;

% suspension 8<->9 + tire
K9(8,8) = K9(8,8) + k_s;   K9(8,9) = K9(8,9) - k_s;
K9(9,8) = K9(9,8) - k_s;   K9(9,9) = K9(9,9) + k_s + k_t;

C9 = zeros(9);
C9(8,8)=C9(8,8)+c_s; C9(8,9)=C9(8,9)-c_s;
C9(9,8)=C9(9,8)-c_s; C9(9,9)=C9(9,9)+c_s+c_t;


%% 3. ISO 8608 CLASS B RANDOM ROAD PROFILE


Fs = 1000;
T = 5;
t_input = 0:1/Fs:T-1/Fs;
v = 20/3.6;
L = v*T*1.05;
xline = linspace(0,L,20000);

Gd0 = 64e-6; n0 = 0.1; w = 2;
dn = 1/L;
n = (1:4000)*dn;
Gd = Gd0*(n/n0).^(-w);
phi = 2*pi*rand(size(n));

h = sum(sqrt(2*dn*Gd(:)).*cos(2*pi*(n(:)*xline)+phi(:)),1);

y = interp1(xline,h,v*t_input,'linear','extrap');
ydot = gradient(y,1/Fs);

F = zeros(9,length(t_input));
F(9,:) = k_t*y + c_t*ydot;

% y = interp1(xline,h,v*t_input,'linear','extrap');
% ydot = gradient(y,1/Fs);
% F(9,:) = k_t*y + c_t*ydot;

%% PLOT ROAD DISPLACEMENT y(t) AND FORCE F(t)

figure;
plot(t_input, y, 'b', 'LineWidth', 1.0);
xlabel('Time (s)');
ylabel('Road displacement y(t) [m]');
title('ISO 8608 Class B Road Profile - Time Domain');
grid on;

figure;
plot(t_input, F(9,:), 'r', 'LineWidth', 1.0);
xlabel('Time (s)');
ylabel('Excitation force F(t) [N]');
title('Tire Road Force Input F(t) = k_t*y + c_t*y˙');
grid on;


%% 4. INMAN MASS NORMALIZATION & MODAL SOLUTION

L = chol(M9,'lower');
invL = inv(L);

K_hat = invL*K9*invL';
C_hat = invL*C9*invL';
F_hat = invL*F;

[U,D] = eig(K_hat);
omega = sqrt(diag(D));
[omega,idx] = sort(omega);
U = U(:,idx);

Phi_full = inv(L')*U;

% mass-normalize eigenvectors
for i = 1.size(Phi_full,2)
    Phi_full(:,i) = Phi_full(:,i) / sqrt(Phi_full(:,i)'*M9*Phi_full(:,i));
end

fprintf("\nNatural frequencies (Hz):\n");
f = omega/(2*pi);          % Hz
format long g              % nicer numeric formatting
disp(f);


% mode reduction to r DOF
r = 4;
Phi_r = Phi_full(:,1:r);
Omega_r = diag(omega(1:r));

% modal forces using physical F (NOT F_hat)
F_modal = Phi_r' * F;

% modal damping matrix C_m = ΦᵀCΦ
C_modal = Phi_r' * C9 * Phi_r;

% modal damping ratios
zeta = diag(C_modal)./(2*diag(Omega_r));
fprintf("\nModal zeta:\n"); disp(zeta);

%% Duhamel integration
Nt = length(t_input);
dt = t_input(2)-t_input(1);

q = zeros(r,Nt);

for i = 1:r
    wi  = Omega_r(i,i);
    zet = zeta(i);
    wdi = wi*sqrt(1-zet^2);
    Qi  = F_modal(i,:);

    for k = 2:Nt
        tau = 0:dt:t_input(k);
        kernel = exp(-zet*wi*(t_input(k)-tau)) .* sin(wdi*(t_input(k)-tau));
        q(i,k) = sum(Qi(1:k).*kernel) * dt / wdi;
    end
end

x_rec = (Phi_r*q).';


%% 5. PLOTS
figure;
plot(t_input,1e3*x_rec(:,7),'LineWidth',1.2); hold on;
plot(t_input,1e3*x_rec(:,8),'--','LineWidth',1.2);
plot(t_input,1e3*x_rec(:,9),':','LineWidth',1.2);
xlabel('Time (s)'); ylabel('Displacement (mm)');
legend('PEMFC base','Sprung mass','Unsprung mass'); grid on;
title('Response with ISO 8608 Class B Road Input (Inman modal reduction)');


%%  6. MODE SHAPE PLOTS (FULL 9-DOF SYSTEM)
% Normalize each mode for plotting (max abs = 1 per column)
Phi9_plot = Phi_full ./ max(abs(Phi_full), [], 1);

% Optional: DOF labels for nicer x-ticks
dof_labels = { ...
    'EP top', ...      % 1
    'BPP+',   ...      % 2
    'GDL+',   ...      % 3
    'PEM',    ...      % 4
    'GDL-',   ...      % 5
    'BPP-',   ...      % 6
    'EP bottom', ...   % 7
    'Sprung mass', ... % 8
    'Unsprung mass'};  % 9

n_modes_to_plot = 6;   % plot first 6 modes (change if you want more)

figure;
for k = 1:n_modes_to_plot
    subplot(3,2,k);
    plot(1:9, Phi9_plot(:,k), '-o', 'LineWidth', 1.5);
    grid on;
    xlim([1 9]);
    xticks(1:9);
    xticklabels(dof_labels);
    xtickangle(45);
    ylabel('Normalized amplitude');
    title(sprintf('Mode %d, f = %.2f Hz', k, omega(k)/(2*pi)));
end
sgtitle('Full 9-DOF PEMFC + Suspension Mode Shapes');
