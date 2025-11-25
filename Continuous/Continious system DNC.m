%% ============================================================
%   Continuous 7-Layer PEMFC PDE Model (Explicit Central Difference)
%
%   Governing PDE:
%       ρ(x) A(x) w_tt = ( EA(x) w_x )_x
%
%   Boundary Conditions (simplified as springs + dampers + road input):
%       x = 0 :  F = -k_top * w(0) - c_top * w_t(0)
%       x = L :  F =  k_bot * (y - w(L)) + c_bot * (y_t - w_t(L))
%
%   This script includes:
%       1. Seven-layer material configuration
%       2. Spatial discretization (finite difference)
%       3. Explicit central-difference time integration
%       4. Base excitation y(t) applied at the bottom boundary
% ============================================================

clear; clc; close all;

%% ------------------------------------------------------------
% 1. Seven-layer stack parameters (thickness / E / density / area)
%% ------------------------------------------------------------
t_layer  = [0.010; 0.002; 0.0008; 0.0002; 0.0008; 0.002; 0.010];
E_layer  = [2e11; 5e10; 1.4286e9; 2e8; 1.4286e9; 5e10; 2e11];
rho_layer = [7800; 4500; 1575; 2000; 1575; 4500; 7800];
A_layer   = 0.01 * ones(7,1);

%% ------------------------------------------------------------
% 2. Boundary parameters (simplified suspension)
%% ------------------------------------------------------------
% Top (free boundary)
k_top = 0;
c_top = 0;

% Bottom: total stiffness/damping of seat + suspension + tire
k_bot = 80000 + 500000;    % 580000 N/m
c_bot = 350 + 15020;       % 15370 Ns/m

%% ------------------------------------------------------------
% 3. Spatial discretization
%% ------------------------------------------------------------
L  = sum(t_layer);
Nx = 301;
dx = L / (Nx - 1);
x  = linspace(0, L, Nx).';

E    = zeros(Nx,1);
rhoA = zeros(Nx,1);
EA   = zeros(Nx,1);

x0 = 0;
for k = 1:7
    x1  = x0 + t_layer(k);
    idx = (x >= x0 & x <= x1);

    E(idx)    = E_layer(k);
    rhoA(idx) = rho_layer(k) * A_layer(k);
    EA(idx)   = E_layer(k) * A_layer(k);

    x0 = x1;
end

EA_edge = 0.5 * (EA(1:end-1) + EA(2:end));

%% ------------------------------------------------------------
% 4. CFL time-step size (stability limit)
%% ------------------------------------------------------------
c_max = max( sqrt(E_layer ./ rho_layer) );
dt     = 0.20 * dx / c_max;

T_end = 0.02;
Nt    = floor(T_end / dt);

fprintf("\n dx = %.3e,   dt = %.3e,   Nt = %d\n", dx, dt, Nt);

%% ------------------------------------------------------------
% 5. Base excitation y(t)
%% ------------------------------------------------------------
A_road = 0.001;
f_exc  = 200;
w_exc  = 2*pi*f_exc;

y_fun  = @(t) A_road * sin(w_exc*t);
yd_fun = @(t) A_road * w_exc * cos(w_exc*t);

%% ------------------------------------------------------------
% 6. Initial conditions
%% ------------------------------------------------------------
w_prev = zeros(Nx,1);
w_curr = zeros(Nx,1);
v_curr = zeros(Nx,1);

%% ------------------------------------------------------------
% 7. Storage for results
%% ------------------------------------------------------------
save_every = max(1, floor(Nt/200));
nsave = floor(Nt / save_every) + 1;

w_hist = zeros(Nx, nsave);
t_hist = zeros(1,  nsave);

idx_save = 1;
w_hist(:,idx_save) = w_curr;
t_hist(idx_save)   = 0;

%% ------------------------------------------------------------
% 8. Time integration loop (explicit central difference)
%% ------------------------------------------------------------
for n = 1:Nt

    t_n = n * dt;

    RHS_n = compute_RHS(w_curr, v_curr, EA_edge, dx, ...
                        k_top, c_top, k_bot, c_bot, ...
                        y_fun(t_n), yd_fun(t_n));

    % Explicit central-difference update
    w_next = 2*w_curr - w_prev + dt^2 * (RHS_n ./ rhoA);
    v_next = (w_next - w_curr) / dt;

    % Shift states
    w_prev = w_curr;
    w_curr = w_next;
    v_curr = v_next;

    % Save data
    if mod(n, save_every) == 0
        idx_save = idx_save + 1;
        w_hist(:,idx_save) = w_curr;
        t_hist(idx_save)   = t_n;
    end
end

%% ------------------------------------------------------------
% 9. Final displacement shape
%% ------------------------------------------------------------
figure; 
plot(x, w_curr, 'LineWidth', 2);
xlabel('x (m)');
ylabel('w (m)');
title('Final displacement distribution w(x)');
grid on;

%% ============================================================
% Local function: RHS builder
% ============================================================
function RHS = compute_RHS(w, v, EA_edge, dx, ...
                           k_top, c_top, k_bot, c_bot, ...
                           y, yd)

Nx  = length(w);
RHS = zeros(Nx,1);

% Internal nodes: (EA w_x)_x
for i = 2:Nx-1
    F_L = EA_edge(i-1)*(w(i)-w(i-1))/dx;
    F_R = EA_edge(i)  *(w(i+1)-w(i))/dx;
    RHS(i) = (F_R - F_L) / dx;
end

% Top boundary x = 0
F_left  = 0;
F_right = EA_edge(1)*(w(2)-w(1))/dx;
RHS(1)  = (F_right - F_left)/dx;

% Bottom boundary x = L
F_leftN  = EA_edge(end)*(w(end)-w(end-1))/dx;
F_rightN = k_bot*(y - w(end)) + c_bot*(yd - v(end));
RHS(end) = (F_rightN - F_leftN)/dx;

end


%% ============================================================
%   PEMFC Continuous Model - Axial Mode Shapes
%   Compute natural frequencies + mode shapes φ_n(x)
%   (Using finite difference: K φ = ω² M φ)
% ============================================================

%% --------------------------
%  1. Geometry (your 7 layers)
% --------------------------
t_layer = [0.010; 0.002; 0.0008; 0.0002; 0.0008; 0.002; 0.010];
E_layer = [2e11; 5e10; 1.4286e9; 2e8; 1.4286e9; 5e10; 2e11];
rho_layer = [7800; 4500; 1575; 2000; 1575; 4500; 7800];
A_layer = 0.01 * ones(7,1);

L = sum(t_layer);
Nx = 300;                % Number of grid points
dx = L/(Nx-1);
x = linspace(0,L,Nx).';

%% --------------------------
% 2. Build EA(x), rhoA(x)
% --------------------------
EA   = zeros(Nx,1);
rhoA = zeros(Nx,1);

x0 = 0;
for k = 1:length(t_layer)
    x1 = x0 + t_layer(k);
    idx = (x >= x0 & x <= x1);

    EA(idx)   = E_layer(k) * A_layer(k);
    rhoA(idx) = rho_layer(k) * A_layer(k);

    x0 = x1;
end

%% --------------------------
% 3. Assemble stiffness K, mass M
% --------------------------
M = diag(rhoA * dx);
K = zeros(Nx, Nx);

for i = 2:Nx-1
    EA_L = (EA(i) + EA(i-1))/2;
    EA_R = (EA(i) + EA(i+1))/2;

    K(i,i-1) = -EA_L/dx;
    K(i,i)   = (EA_L + EA_R)/dx;
    K(i,i+1) = -EA_R/dx;
end

% Boundary x=0 (spring)
k_top = 3e5;
EA_R = (EA(1) + EA(2))/2;
K(1,1) = EA_R/dx + k_top;
K(1,2) = -EA_R/dx;

% Boundary x=L (spring)
k_bot = 3e5;
EA_L = (EA(Nx) + EA(Nx-1))/2;
K(Nx,Nx-1) = -EA_L/dx;
K(Nx,Nx)   = EA_L/dx + k_bot;

%% --------------------------
% 4. Solve K φ = ω² M φ
% --------------------------
n_modes = 4;
[Phi, Omega2] = eigs(K, M, n_modes, 'smallestabs');
omega_n = sqrt(diag(Omega2));
f_n = omega_n/(2*pi);

% Normalize
for k = 1:n_modes
    Phi(:,k) = Phi(:,k) / max(abs(Phi(:,k)));
end

%% --------------------------
% 5. Plot mode shapes
% --------------------------
figure; hold on;
colors = lines(n_modes);

for k = 1:n_modes
    plot(x, Phi(:,k), 'LineWidth', 2, 'Color', colors(k,:));
end

xlabel('x position (m)');
ylabel('\phi_n(x)  (normalized)');
title('PEMFC Continuous Model - Axial Mode Shapes');

legend( ...
    ['Mode 1 (' num2str(f_n(1),'%.1f') ' Hz)'], ...
    ['Mode 2 (' num2str(f_n(2),'%.1f') ' Hz)'], ...
    ['Mode 3 (' num2str(f_n(3),'%.1f') ' Hz)'], ...
    ['Mode 4 (' num2str(f_n(4),'%.1f') ' Hz)'], ...
    'Location','Best' ...
);

grid on;

