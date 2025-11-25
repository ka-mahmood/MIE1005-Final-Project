%% ==============================================
%   1. System Parameters (BUS data fully substituted)
% ==============================================
% Masses [kg]
m_end  = 100;       % Upper plate (driver–seat mass)
m_cell = 0.10;      % PEMFC core mass
m_s    = 4500;      % Sprung mass (bus body on front axle)
m_us   = 500;       % Unsprung mass (wheel–axle assembly)

% Stiffness [N/m]
k_clamp = 25000;    % Driver–seat spring stiffness
k_cell  = 8.0e7;    % Equivalent stiffness of PEMFC core
k_end   = 300000;   % Front axle bellows stiffness
k_s     = 300000;   % Suspension spring stiffness
k_t     = 1600000;  % Tire stiffness

% Damping [N·s/m]
c_s  = 1000;        % Suspension damper (seat–body)
c_us = 20000;       % Front axle damper
c_t  = 150;         % Tire damping

% Degree of Freedom (DOF) ordering:
% x1 = Upper plate (driver–seat)
% x2 = PEMFC core
% x3 = Lower plate (bus base)
% x4 = Sprung mass (vehicle body)
% x5 = Unsprung mass (wheel/axle)
ndof = 5;

%% ==============================================
%   2. Mass Matrix
% ==============================================
M = diag([m_end, m_cell, m_end, m_s, m_us]);

%% ==============================================
%   3. Stiffness Matrix  (corrected with k_clamp coupling)
% ==============================================
K = [ k_clamp + k_cell,   -k_cell,           -k_clamp,       0,          0;
      -k_cell,            2*k_cell,          -k_cell,        0,          0;
      -k_clamp,           -k_cell,   k_clamp + k_cell + k_end, -k_end,    0;
       0,                  0,                -k_end,     k_end + k_s,    -k_s;
       0,                  0,                 0,             -k_s,    k_s + k_t ];

%% ==============================================
%   4. Damping Matrix (suspension and tire branch)
% ==============================================
C = zeros(ndof);
C(4,4) =  c_s;     C(4,5) = -c_s;
C(5,4) = -c_s;     C(5,5) =  c_s + c_us;

%% ==============================================
%   5. Road Excitation (ISO 8608 Class B)
% ==============================================
Fs = 1000;              % Sampling frequency [Hz]
T = 5;                  % Simulation duration [s]
t = 0:1/Fs:T-1/Fs;      % Time vector
N = length(t);          % Number of time steps

v = 20/3.6;             % Vehicle speed [m/s]
L = v*T*1.05;           % Road length (5% margin)
xline = linspace(0, L, 20000);   % Spatial sampling points

% ISO 8608 parameters (Class B)
Gd0 = 64e-6;   % Reference PSD [m^3]
n0  = 0.1;     % Reference spatial frequency [1/m]
w   = 2;       % Waviness exponent

% Spatial frequency grid
dn = 1/L; 
n  = (1:4000)*dn;                   % Spatial frequency [1/m]
Gd = Gd0 * (n/n0).^(-w);            % Spatial PSD [m^3]
phi = 2*pi*rand(size(n));           % Random phase

% Synthesize spatial road profile h(x)
h = sum( sqrt(2*dn*Gd(:)) .* cos(2*pi*(n(:)*xline) + phi(:)), 1 );

%% --- Spatial → Temporal Mapping: y(t) = h(v*t) ---
y  = interp1(xline, h, v*t, 'linear', 'extrap');  % Road displacement input [m]
yd = gradient(y, 1/Fs);                           % Road velocity input [m/s]

%% --- External Force Vector (acts only on DOF 5: wheel/axle) ---
f = zeros(ndof, N);
f(5,:) = k_t * y + c_us * yd;

%% ==============================================
%   6. Time Integration (Newmark–β, average acceleration method)
% ==============================================
beta  = 1/4;  
gamma = 1/2;    % Unconditionally stable parameters

x   = zeros(ndof, N);   % Displacement
x_d = zeros(ndof, N);   % Velocity
a   = zeros(ndof, N);   % Acceleration

% Initial acceleration
a(:,1) = M \ (f(:,1) - C*x_d(:,1) - K*x(:,1));

% Effective stiffness matrix
dt = 1/Fs;
K_eff = K + (gamma/(beta*dt))*C + (1/(beta*dt^2))*M;
invK_eff = inv(K_eff);

% --- Time-stepping loop ---
for i = 1:N-1

    % Predictor step
    x_pred = x(:,i) + dt*x_d(:,i) + (0.5 - beta)*dt^2 * a(:,i);
    v_pred = x_d(:,i) + (1 - gamma)*dt * a(:,i);

    % Effective force
    f_eff = f(:,i+1) ...
            + M*((1/(beta*dt^2))*x_pred) ...
            + C*((gamma/(beta*dt))*x_pred);

    % Corrector step
    x(:,i+1) = invK_eff * f_eff;
    a(:,i+1) = (1/(beta*dt^2)) * (x(:,i+1) - x_pred);
    x_d(:,i+1) = v_pred + gamma*dt * a(:,i+1);
end

%% ==============================================
%   7. Internal Force Computation
% ==============================================
F_clamp  = k_clamp * x(1,:);                  % Clamp force
F_core   = k_cell * (x(1,:) - x(2,:));        % PEMFC core compression
F_end    = k_end  * (x(3,:) - x(4,:));        % Lower plate–vehicle force
F_spring = k_s    * (x(4,:) - x(5,:));        % Suspension spring force
F_damper = c_s    * (x_d(4,:) - x_d(5,:));    % Suspension damper force
F_tire   = k_t * (x(5,:) - y) + c_us * (x_d(5,:) - yd);  % Tire–road force

%% ==============================================
%   8. Results – Time Responses
% ==============================================
figure;
plot(t, x(4,:)*1000, 'b', 'LineWidth',1.2); hold on;
plot(t, x(5,:)*1000, 'k', 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('Displacement (mm)');
legend('Sprung mass','Unsprung mass');
title('Vehicle Response with Detailed PEMFC Stack (5-DOF)');
grid on;

% --- Time response summary (Displacement, Velocity, Acceleration) ---
figure;
subplot(3,1,1);
plot(t, x(4,:)*1000, 'b', 'LineWidth',1.2); hold on;
plot(t, x(5,:)*1000, 'k', 'LineWidth',1.2);
ylabel('Displacement (mm)');
legend('Sprung mass','Unsprung mass');
title('Time Response of Vehicle Suspension');

subplot(3,1,2);
plot(t, x_d(4,:), 'b', 'LineWidth',1.2); hold on;
plot(t, x_d(5,:), 'k', 'LineWidth',1.2);
ylabel('Velocity (m/s)');
legend('Sprung mass','Unsprung mass');

subplot(3,1,3);
plot(t, a(4,:), 'b', 'LineWidth',1.2); hold on;
plot(t, a(5,:), 'k', 'LineWidth',1.2);
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');
legend('Sprung mass','Unsprung mass');
grid on;

%% --- Internal PEMFC Layer Displacement ---
figure;
plot(t, x(1:5,:)'*1000, 'LineWidth',1);
xlabel('Time (s)'); ylabel('Displacement (mm)');
legend('BPP+','GDL+','PEM','GDL−','BPP−');
title('Internal PEMFC Layer Vibrations');
grid on;

%% ==============================================
%   9. Modal Analysis
% ==============================================
[Phi, Omega2] = eig(K, M);          % Solve K*phi = ω²*M*phi
omega_n = sqrt(diag(Omega2));       % Natural angular frequencies [rad/s]
f_n = omega_n / (2*pi);             % Natural frequencies [Hz]

% Mass normalization (Phiᵀ*M*Phi = I)
for k = 1:length(f_n)
    mk = Phi(:,k)' * M * Phi(:,k);
    Phi(:,k) = Phi(:,k) / sqrt(mk);
end

% Print modal table
fprintf('\n================ Modal Analysis ================\n');
fprintf('Mode\tFreq(Hz)\tDominant DOF\n');
dof_names = ["Upper plate","PEMFC core","Lower plate","Sprung mass","Wheel/axle"];
for k = 1:ndof
    [~, idx] = max(abs(Phi(:,k)));
    fprintf('%d\t%.3f\t\t%s\n', k, f_n(k), dof_names(idx));
end

% --- Plot normalized mode shapes (full 5-DOF system) ---
figure;
for k = 1:ndof
    subplot(ceil(ndof/2), 2, k);
    plot(1:ndof, Phi(:,k)/max(abs(Phi(:,k))), '-o', 'LineWidth', 1.2);
    xlabel('Degree of Freedom');
    ylabel('Normalized amplitude');
    title(['Mode ', num2str(k), ', f = ', num2str(f_n(k), '%.2f'), ' Hz']);
    grid on;
end
sgtitle('Full 5-DOF PEMFC + Suspension Mode Shapes');

%% ==============================================
%   10. Natural Frequency Summary
% ==============================================
fprintf('\n=== Natural Frequencies (Hz) ===\n');
for k = 1:length(f_n)
    fprintf('Mode %d:  %.3f Hz\n', k, f_n(k));
end
