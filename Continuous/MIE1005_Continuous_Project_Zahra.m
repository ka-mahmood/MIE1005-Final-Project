clear; close all; clc;

%% 1. PHYSICAL PARAMETERS FOR EACH REAL PART (EDIT THESE)

% Endplates
E_ep   = 190e9;      % Pa
rho_ep = 0.8064/(0.02*5040*(0.001^2));        % kg/m^3
A_ep   = 5040*(0.001^2);                      % m^2
L_ep   = 0.02;                                % m  (thickness of ONE endplate)

% BPP ribs
E_rib   = 105e9;     % Pa
rho_rib = (22.275-0.495)*0.001/((110+550*2)*(0.001^2)*0.001);  % kg/m^3
A_rib   = (110+550*2)*(0.001^2);             % m^2

% BPP main plate (solid plate the ribs sit on)
E_bpp_main   = 105e9;   % Pa
rho_bpp_main = rho_rib; % same material
A_bpp_main   = 1650*(0.001^2);   % m^2 (approx active area)
L_bpp_main   = 0.003;           % m
L_bpp_rib    = 0.001;           % m
L_bpp        = 0.001;           % EFFECTIVE axial thickness of one BPP (as before)

% BPP channels (void or coolant)
E_ch   = 1.0e6;       % Pa (very soft, almost zero)
rho_ch = 0;           % kg/m^3 (fluid) or ~0
A_ch   = 110*1*(0.001^2);       % m^2

% GDL- (negative side)
E_gdl_minus   = 117e9;    % Pa
rho_gdl_minus = 0.002937/(0.0002*1650*(0.001^2));     % kg/m^3
A_gdl         = 1650*(0.001^2);     % m^2
L_gdl         = 0.0002;             % m

% GDL+ (positive side)
E_gdl_plus   = 10e9;     % Pa
rho_gdl_plus = 0.000145/(0.0002*1650*(0.001^2));      % kg/m^3

% Gaskets (same material both sides, two per side)
E_g   = 50e6;           % Pa
rho_g = 0.4185*0.5*0.001/(0.00045*275*(0.001^2)); % kg/m^3
A_g   = 275*(0.001^2);  % m^2 (one gasket)
L_g   = 0.00045;        % m

% PEM
E_pem   = 241e6;        % Pa
rho_pem = 0.000165/(0.00005*1650*(0.001^2));   % kg/m^3
A_pem   = 1650*(0.001^2);                      % m^2
L_pem   = 0.00005;                             % m

% Top spring stiffness (clamping frame to vehicle)
k_mount_top = 1e6;
k_endplate  = 4.788e10;
k_top       = (k_endplate * k_mount_top)/(k_endplate + k_mount_top);

% Number of modes to keep in PEMFC continuous model
Nmodes = 4;

%%% CHANGED: no internal PEMFC damping; we will NOT use zeta here anymore

% Suspension / quarter-car parameters (same spirit as  9-DOF model)
m_s  = 2500;
m_us = 200;

k_s  = 80000;
c_s  = 350;
k_t  = 0.5e6;
c_t  = 15020;

%% 2. BUILD 7 REGIONS: EA_i, rhoA_i, L_i  (exact PEMFC layout)

% Region 1: bottom endplate
EA1   = E_ep   * A_ep;
rhoA1 = rho_ep * A_ep;
L1    = L_ep;

% Region 2: bottom BPP (ribs + channels in PARALLEL axially)
EA2   = E_rib*A_rib + E_ch*A_ch;
rhoA2 = rho_rib*A_rib + rho_ch*A_ch;
L2    = L_bpp;

% Region 3: GDL- + TWO bottom gaskets (parallel)
EA3   = E_gdl_minus*A_gdl + 2*E_g*A_g;
rhoA3 = rho_gdl_minus*A_gdl + 2*rho_g*A_g;
L3    = L_gdl + 2*L_g;

% Region 4: PEM only
EA4   = E_pem*A_pem;
rhoA4 = rho_pem*A_pem;
L4    = L_pem;

% Region 5: GDL+ + TWO top gaskets (parallel)
EA5   = E_gdl_plus*A_gdl + 2*E_g*A_g;
rhoA5 = rho_gdl_plus*A_gdl + 2*rho_g*A_g;
L5    = L_gdl + 2*L_g;

% Region 6: top BPP
EA6   = EA2;
rhoA6 = rhoA2;
L6    = L_bpp;

% Region 7: top endplate
EA7   = EA1;
rhoA7 = rhoA1;
L7    = L_ep;

EA_regions   = [EA1 EA2 EA3 EA4 EA5 EA6 EA7];
rhoA_regions = [rhoA1 rhoA2 rhoA3 rhoA4 rhoA5 rhoA6 rhoA7];
L_regions    = [L1  L2  L3  L4  L5  L6  L7];

Nr = numel(EA_regions);


%% 3. 1D BAR FEM DISCRETIZATION (linear elements)
nel_region = [4 4 4 4 4 4 4];  % elements per region

Ne = sum(nel_region);      % total elements
Nn = Ne + 1;               % total nodes

x      = zeros(Nn,1);
EA_e   = zeros(Ne,1);
rhoA_e = zeros(Ne,1);

e_global = 0;
x_curr   = 0;

for r = 1:Nr
    n_el_r = nel_region(r);
    Lr     = L_regions(r);
    Le     = Lr / n_el_r;
    
    for e = 1:n_el_r
        e_global = e_global + 1;
        
        i = e_global;       % left node index
        j = e_global + 1;   % right node index
        
        if e == 1 && r == 1
            x(i) = x_curr;
        end
        x(j) = x(i) + Le;
        
        EA_e(e_global)   = EA_regions(r);
        rhoA_e(e_global) = rhoA_regions(r);
    end
    
    x_curr = x(e_global+1);
end

L_total = x(end);

% Assemble global K and M
K = zeros(Nn);
M = zeros(Nn);

for e = 1:Ne
    i = e;
    j = e + 1;
    
    Le   = x(j) - x(i);
    EA   = EA_e(e);
    rhoA = rhoA_e(e);
    
    ke = (EA/Le) * [ 1 -1; -1  1];
    me = (rhoA*Le/6) * [2 1; 1 2];
    
    dofs = [i j];
    K(dofs,dofs) = K(dofs,dofs) + ke;
    M(dofs,dofs) = M(dofs,dofs) + me;
end

% Boundary conditions for eigenproblem:
% - Bottom x=0: free  -> no modification
% - Top x=L: spring k_top in parallel with bar -> add to last DOF
K(end,end) = K(end,end) + k_top;


%% 4. SOLVE EIGENPROBLEM  [K]phi = omega^2 [M]phi
[Phi, D] = eig(K,M);
omega_all = sqrt(diag(D));
[omega_all, idx] = sort(omega_all);
Phi = Phi(:,idx);

omega_n = omega_all(1:Nmodes);
fn      = omega_n/(2*pi);

disp('Natural frequencies of layered PEMFC stack (Hz):');
disp(fn.');

% Mass-normalize modes: Phi' * M * Phi = I
for n = 1:Nmodes
    mn = sqrt(Phi(:,n)' * M * Phi(:,n));
    Phi(:,n) = Phi(:,n) / mn;
end
%% MODE SHAPE PLOTS
figure;
hold on;

colors = lines(Nmodes);   % color set

for n = 1:Nmodes
    plot(x, Phi(:,n), 'LineWidth', 2, 'Color', colors(n,:));
end

grid on;
xlabel('x position (m)');
ylabel('\phi_n(x)  (normalized)');
title('PEMFC Continuous Model - Axial Mode Shapes');
legend(arrayfun(@(k) sprintf('Mode %d (%.1f Hz)', k, fn(k)), 1:Nmodes, 'UniformOutput', false));


%% 5. ROAD PROFILE y(t)  (ISO 8608 Class B) AND SUSPENSION RESPONSE
Fs = 1000;      % sampling frequency [Hz]
T  = 10;         % duration [s]
t  = 0:1/Fs:T-1/Fs;
Nt = numel(t);
dt = 1/Fs;

v     = 20/3.6;          % vehicle speed [m/s]
Lroad = v*T*1.05;        % road length
xline = linspace(0,Lroad,20000);

% ISO 8608 Class B spectral parameters
Gd0 = 16e-6;
n0  = 0.1;
w   = 2;
dn  = 1/Lroad;
n   = (1:4000)*dn;
Gd  = Gd0*(n/n0).^(-w);
phi = 2*pi*rand(size(n));

% Spatial road profile h(x)
h = sum( sqrt(2*dn*Gd(:)).*cos(2*pi*(n(:)*xline) + phi(:)), 1 );

% Map to time domain y(t)
y    = interp1(xline, h, v*t, 'linear', 'extrap');  % road displacement
ydot = gradient(y, dt);                             % road velocity

% State variables: [x_s; x_us; v_s; v_us]
xs  = 0;   xus  = 0;
vs  = 0;   vus  = 0;

F_eq = zeros(1,Nt);

for k = 1:Nt-1
    yk    = y(k);
    ydotk = ydot(k);
    
    % accelerations from quarter-car equations
    a_s  = (-c_s*(vs - vus) - k_s*(xs - xus)) / m_s;
    a_us = (-c_s*(vus - vs) - k_s*(xus - xs) ...
            - c_t*(vus - ydotk) - k_t*(xus - yk)) / m_us;
    
    % semi-implicit Euler integration (stable enough for dt=1e-3)
    vs  = vs  + a_s*dt;
    vus = vus + a_us*dt;
    xs  = xs  + vs*dt;
    xus = xus + vus*dt;
    
    % Suspension force transmitted to sprung mass (and hence to PEMFC mount)
    F_eq(k) = k_s*(xs - xus) + c_s*(vs - vus);
end
F_eq(Nt) = F_eq(Nt-1);   % last sample

% NOTE: F_eq(t) now depends on y(t) *and* includes c_s, c_t effects.
% No damping will be added inside the PEMFC model itself.


%% 6. MODAL RESPONSE (UNDAMPED Inman Duhamel integral 6.8)
%    q¨_n + omega_n^2 q_n = Phi(1,n) * F_eq(t) -- WINDOW 4.8

q = zeros(Nmodes,Nt);

for n = 1:Nmodes
    w  = omega_n(n);
    Gn = Phi(1,n);          % generalized force coefficient at bottom DOF
    
    for k = 2:Nt
        tau    = t(1:k);
        kernel = sin(w*(t(k)-tau));     % UNDAMPED kernel
        q(n,k) = Gn/w * trapz(tau, F_eq(1:k) .* kernel);
    end
end

% Reconstruct physical displacement: u(x,t) = sum_n Phi_n(x) q_n(t)
u = (Phi(:,1:Nmodes) * q).';   % size Nt x Nn


%% 7. PLOTS
% Bottom and top endplate responses
figure;
plot(t, 1e3*u(:,1), 'LineWidth', 1.5); hold on;
plot(t, 1e3*u(:,4),'LineWidth', 1.5);
plot(t, 1e3*u(:,end), '--', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Displacement (mm)');
legend('Bottom endplate (x=0)','Top endplate (x=L)');
title('Layered PEMFC axial response (continuous, undamped PEMFC; damped suspension)');

% Full x-t field
figure;
surf(x, t, 1e3*u, 'EdgeColor','none');
xlabel('x (m)');
ylabel('t (s)');
zlabel('u(x,t) [mm]');
title('Continuous layered PEMFC stack vibration (axial)');
colormap turbo;
view(30,30);


%%

u_bottom_m  = u(:,1);           % meters
u_bottom_mm = 1e3*u_bottom_m;   % mm

% Design low-pass filter
fc = 20;                        % Hz
[b,a] = butter(4, fc/(Fs/2));

% Filter in mm units
u_lp_mm = filtfilt(b,a, u_bottom_mm);

figure;
plot(t, u_bottom_mm, 'Color',[0.8 0.4 0.4]); hold on;
plot(t, u_lp_mm, 'k', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Displacement (mm)');
legend('Original (unfiltered)','Low-pass filtered (20 Hz)');
title('Low-pass filtered PEMFC response');

% Optional: turn off scientific notation on the y-axis
ax = gca;
ax.YAxis.Exponent = 0;

