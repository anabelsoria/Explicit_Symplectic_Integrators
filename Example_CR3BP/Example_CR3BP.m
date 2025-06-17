%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Example Script for Symplectic Integrations
% 
% This script compares symplectic (SI) and Runge-Kutta (RK) integrators by
% propagating an orbit in the Circular Restricted Three-Body Problem
% (CR3BP) over several revolutions. It plots the 3D trajectory and the 
% Jacobi constant variation over time to evaluate energy preservation.
%
% Author: Anabel Soria-Carro 
% Date:   May 22, 2025
% Affiliation: The University of Texas at Austin
%              Controls Group for Distributed and Uncertain Systems (CDUS)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

clear; clc; close all;

% Add all subfolders of the parent directory to the path
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')))

%% ====================== Data Setup ======================

orbit_type = 'NRHO_L2_S';  % Specify orbit type (DRO, NRHO_L2_S, Halo_L1_N)
center     = 'bary';
p = CR3BPOrbit(orbit_type, center);

Nrevs = 1;          % Number of revolutions to propagate
Nsteps = 1000;       % Number of time steps per revolution

order  = 6;         % Integrators order
scheme = 2;         % Störmer-Verlet scheme 1 or 2

%% ====================== Propagate Nrevs ======================

% -------------------- SYMPLECTIC INTEGRATOR --------------------
% Create instance of the symplectic integrator class
SI_obj = SI(p, order, scheme);  

% Define propagation parameters
t0 = 0;                      % Initial time
tf = Nrevs * p.Tp;           % Final time = Nrevs full orbital periods
dt = p.Tp / Nsteps;          % Step size

% Propagate using symplectic integrator
tic
[X_SI, t_SI] = SI_obj.propagate(t0, tf, dt);
toc

% ------------------------ RUNGE-KUTTA ------------------------
% Create instance of the Runge-Kutta integrator class
RK_obj = RK(p, order);

% Propagate using RK integrator
tic
[X_RK, t_RK] = RK_obj.propagate(p.nu0, t0, tf, dt, ...
                                @(t, x) p.DS.Hamiltons_EOM(t, x));
toc
%% ======================== POST-PROCESSING =========================

% ------------------------ Plot Orbits -----------------------------
figure; hold on; grid on; axis equal;
plot3(X_SI(1,:), X_SI(2,:), X_SI(3,:), LineWidth=2); 
plot3(X_RK(1,:), X_RK(2,:), X_RK(3,:), LineWidth=2); 
labels3d('LU');

% ------------------- Plot Jacobi Constant Drift -------------------
% Compute the Jacobi constant over time for both integrators
% and plot the absolute difference from the initial value.

% --- RK Integrator ---
X_RK_cart = p.DS.nu2xi(X_RK);
C_RK       = p.DS.jacobiconstant(X_RK_cart);        
diff_C_RK  = abs(C_RK - C_RK(1));           % Deviation from initial value

% --- Symplectic Integrator ---
X_SI_cart = p.DS.nu2xi(X_SI);
C_SI       = p.DS.jacobiconstant(X_SI_cart);        
diff_C_SI  = abs(C_SI - C_SI(1));           % Deviation from initial value

% --- Plotting the Drift ---
figure; hold on;
plot(t_SI/p.Tp, diff_C_SI, LineWidth=2, DisplayName='SI');    
plot(t_RK/p.Tp, diff_C_RK, LineWidth=2, DisplayName='RK');    
grid on;
xlabel('Revolutions', 'Interpreter', 'latex');                
ylabel('$|C - C_0|$', 'Interpreter', 'latex');                
set(gca, 'yscale', 'log');                                    
legend('Interpreter', 'latex', 'Location','best');                               
