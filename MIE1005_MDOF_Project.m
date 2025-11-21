%% Material properties 

% Endplate
E_e = 190e9;
A_e = 5040e-6;
L_e = 20e-3;
p_e = 8;

% BPP main
E_bpp = 105e9;
A_bpp_main = 1650e-6;
L_bpp_main = 3e-3;

% BPP edge
n_edges = 2;
A_bpp_edge = 550-6*n_edges;
L_bpp_edge = 1e-3;

% BPP rib
n_ribs = 1;
A_bpp_rib = 110e-6*n_ribs;
L_bpp_rib = 1e-3;
p_bpp = 4.5;

% Gasket
E_g = 50e6;
A_g = 275e-6;
L_g = 0.45e-3;
p_g = 1.2;

% GDL - anode *** UPDATE THE SOLIDWORKS TO REFLECT THIS
E_gdl_a = 117e9;
A_gdl_a = A_bpp_main - 2*A_g; 
L_gdl_a = 0.2e-3;
p_gdl_a = 8.9;

% GDL - cathode
E_gdl_c = 10e9;
A_gdl_c = A_bpp_main - 2*A_g;
L_gdl_c = 0.2e-3;
p_gdl_c = 0.44;

% PEM
E_pem = 241e6;
A_pem = A_bpp_main; % 1650e-6; % Need to change to BPP main area for now
L_pem = 0.05e-3;
p_pem = 2.3; % choose higher range for the PEM density

% Bolt properties
E_bolt = E_e;
n_bolt = 4;
r_bolt = 8e-3;
A_bolt = pi*(r_bolt)^2*n_bolt; % Calculate the bolt cross section
L_stack = 2*(L_bpp_edge + L_bpp_main) + L_pem + L_g;

% Properties of bus suspension system
m_s = 2500;
m_us = 350;
k_s = 80000;
k_us = 500000;
b_s = 350;
b_us = 15020;

%% Calculate the spring stiffness for each element
k_e = E_e*A_e/L_e; % Endplate stiffness

% Calculate the summed k_bpp from the main, rib, and edge
k_bpp_main = E_bpp*A_bpp_main/L_bpp_main;
k_bpp_rib = E_bpp*A_bpp_rib/L_bpp_rib;
k_bpp_edge = E_bpp * A_bpp_edge / L_bpp_edge;
k_bpp = (1/(k_bpp_rib + k_bpp_edge) + 1/k_bpp_main)^-1;

% Calculate the other component stiffnesses
k_gdl_a = E_gdl_a * A_gdl_a / L_gdl_a;
k_gdl_c = E_gdl_c * A_gdl_c / L_gdl_c;
k_gasket = E_g * A_g / L_g;
k_pem = E_pem * A_pem / L_pem;
k_clamp = E_bolt*A_bolt / L_stack;

% Treat the subsequent masses as a lump mass
% Calculate the cell stack equivalent stiffness
k_pem_gasket = (2/k_gasket + 1/k_pem)^-1;
k_pem_gdl = (1/k_gdl_a + 1/k_pem + 1/k_gdl_c)^-1;
k_pem_gdl_gasket = k_pem_gasket + k_pem_gdl;

% Add in the BPP to include into the system
k_cell = (1/k_bpp + 1/k_pem_gdl_gasket + 1/k_bpp)^-1;

% Since we are interested in the movement of k_cell, not
%   including the endplate, we will not further collapse the system

%% Calculate the system using the derived values

% Calculate the masses
m_e = L_e*A_e*p_e;
m_bpp = L_bpp_main * A_bpp_main * p_bpp; 
m_g = L_g*A_g*p_g;
m_gdl_a = A_gdl_a*L_gdl_a*p_gdl_a;
m_gdl_c = A_gdl_c*L_gdl_c*p_gdl_c;
m_pem = A_pem*L_pem*p_pem;

% Sum the masses for m_cell
m_cell = m_bpp + m_g + m_gdl_a + m_pem + m_gdl_c + m_bpp;

M = diag([m_us, m_s, m_e, m_cell]);
K = [k_s + k_us, -k_s, 0, 0;
    -k_s, k_s + k_e, -k_e, 0;
    0,  -k_e, k_e + k_cell + k_clamp, -k_cell;
    0, 0, -k_cell, k_cell + k_e];
C = [b_us + b_s, -b_s, 0, 0;
     -b_s, b_s, 0, 0;
     0, 0, 0, 0;
     0, 0, 0, 0];

%% Solve for the mode shapes, natural freq., and damped natural freq.
% Natural modes only apply to undamped system
K_equiv = M^(-1/2)*K*M^(-1/2);
C_equiv = M^(-1/2)*C*M^(-1/2);

% Compute modes and natural frequencies
[V, D] = eig(K_equiv);
w = sqrt(diag(D));
u1 = V(:, 1)/norm(V(:, 1));
u2 = V(:, 2)/norm(V(:, 2));
u3 = V(:, 3)/norm(V(:, 3));
u4 = V(:, 4)/norm(V(:, 4));

% Compute relevant equivalent matrices for stiffness & P
P = horzcat(u1, u2, u3, u4);
KP = P.' * K_equiv * P;
PtP = P.' * P;
CP = P.' * C_equiv * P;

% Compute the damping
d1 =  CP(1,1)/(2*w(1));
d2 =  CP(2,2)/(2*w(2));
d3 =  CP(3,3)/(2*w(3));
d4 =  CP(4,4)/(2*w(4));

% Compute the damped natural frequencies
wd1 =  w(1)*sqrt(1-d1^2);
wd2 =  w(2)*sqrt(1-d2^2);
wd3 =  w(3)*sqrt(1-d3^2);
wd4 =  w(4)*sqrt(1-d4^2);

%% Solve for the displacement of the system given a force

% Solve for the modal equations
% Let the force vector be the sum of the forces from k_us and b_us
syms t
y = @(t) 0.1*sin(10*t); 
ydot = diff(y, t);
F = @(t) ...
    [k_us*y(t) + b_us*ydot(t), 0, 0, 0].';

% Let the initial conditions be 0
% r.. + diag(2*di*wni)*r. + delta_r = PT*M^(-1/2)*F


% Compute the modal force vector
modal_f_mult = P.' * M^(-1/2);
modal_f = modal_f_mult * [1, 0, 0, 0].';
modal_f_forced1 = @(t) modal_f .* F(t);
modal_f_forced2 = @(t) modal_f .* F(t);
modal_f_forced3 = @(t) modal_f .* F(t);
modal_f_forced4 = @(t) modal_f .* F(t);

syms r1(t) r2(t) r3(t) r4(t)

dr1 = diff(r1,t);
dr2 = diff(r2,t);
dr3 = diff(r3,t);
dr4 = diff(r4,t);

ddr1 = diff(dr1,t);
ddr2 = diff(dr2,t);
ddr3 = diff(dr3,t);
ddr4 = diff(dr4,t);
f_x1 = modal_f_forced(1);

eqn1 = ddr1 + diag(2*d1*w(1))*dr1 + w(1)^2*r1 == f_x1(1);
eqn2 = ddr2 + diag(2*d2*w(2))*dr2 + w(2)^2*r2 == 0;
eqn3 = ddr3 + diag(2*d3*w(3))*dr3 + w(3)^2*r3 == 0;
eqn4 = ddr4 + diag(2*d4*w(4))*dr4 + w(4)^2*r4 == 0;

% Solve each ODE equation
cond1 = [r1(0) == 0, dr1(0) == 0];
cond2 = [r2(0) == 0, dr2(0) == 0];
cond3 = [r3(0) == 0, dr3(0) == 0];
cond4 = [r4(0) == 0, dr4(0) == 0];

% Solve for x by converting back from r
convert = M^(-1/2)*P;
x1(t) = convert(:, 1).*dsolve(eqn1, cond1);
x2(t) = convert(:, 2).*dsolve(eqn2, cond2);
x3(t) = convert(:, 3).*dsolve(eqn3, cond3);
x4(t) = convert(:, 4).*dsolve(eqn4, cond4);

t0 = 1;
fplot(x1, [0, 10])
hold on
fplot(x2, [0.5, 10])
fplot(x3, [0.5, 10])
fplot(x4)

%% Plotting
do_plot = true;

if do_plot == true
    % Plot the modes
    figure
    plot(u1)
    hold on
    plot(u2)
    plot(u3)
    plot(u4)
    hold off
    title("Mode Shapes for MDOF PEMFC Model with Suspension")
    ylabel("Normalized Displacement")
    xlabel("Degree of Freedom")
    ax = gca;
    ax.XTick = unique( round(ax.XTick) );
    legend("\omega_1 = " + w(1) + " rad/s", ...
           "\omega_2 = " + w(2) + " rad/s", ...
           "\omega_3 = " + w(3) + " rad/s", ...
           "\omega_4 = " + w(4) + " rad/s", 'Location', 'bestoutside')
end





%% Archive Work
% syms y(t) x1(t) x2(t) x3(t) x4(t)
% ydot = diff(y, t);
% 
% x = [x1, x2, x3, x4];
% xdot = diff(x, t);              
% xddot = diff(xdot, t);
% ode = M*xddot.' + K*x.' + C*xdot.' == 0; % y*[k_us, 0, 0, 0].' + ydot*[b_us, 0, 0, 0].';
% sol = dsolve(ode)
% 
% % Add in the appropriate boundary conditions
% 
% 
% % % Determine the solution numerically
% syms omega
% A = K - omega^2 * M - omega * C;
% eqn = det(A); % natural frequency of the system
% lambda  = vpasolve(eqn == 0);
% eigenval = lambda(lambda>=0); % Use only valid values for natural frequency
% w = sqrt(eigenval);
% 
% 
% % Plot a 2 x 2 set of plots for u1, u2, u3, u4
% % Create a figure for the plots
% figure;
% % Plot the natural frequencies
% plot(u1, 'o-', 'DisplayName', "\omega_1=" + num2str(w_undamp(1)));
% hold on
% plot(u2, 'o-', 'DisplayName', "\omega_2=" + num2str(w_undamp(2)));
% plot(u3, 'o-', 'DisplayName', "\omega_3=" + num2str(w_undamp(3)));
% plot(u4, 'o-', 'DisplayName', "\omega_4=" + num2str(w_undamp(4)));
% hold off
% title('Undamped Natural Frequencies');
% xlabel('Degree of Freedom');
% ylabel('Relative Displacement');
% legend()
% 
% % Plot a 2 x 2 set of plots for u1, u2, u3, u4
% % Create a figure for the plots
% figure;
% % Plot the natural frequencies
% plot(u1, 'o-', 'DisplayName', "\omega_1=" + string(round(w(1), 4)));
% hold on
% plot(u2, 'o-', 'DisplayName', "\omega_2=" + string(round(w(2), 4)));
% plot(u3, 'o-', 'DisplayName', "\omega_3=" + string(round(w(3), 4)));
% plot(u4, 'o-', 'DisplayName', "\omega_4=" + string(round(w(4), 4)));
% hold off
% title('Damped Natural Frequencies  Modal Analysis Not Applicable');
% xlabel('Degree of Freedom');
% ylabel('Relative Displacement');
% legend()