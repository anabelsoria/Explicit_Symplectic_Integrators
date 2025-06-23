%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Example Script for Symplectic Integrations
% 
% This script compares symplectic (SI) and Runge-Kutta (RK) integrators by
% propagating an orbit in the Bicircular Restricted Four-Body Problem (BCR4BP)
% over several revolutions. 
%
% Author: Anabel Soria-Carro 
% Date:   June 20, 2025
% Affiliation: The University of Texas at Austin
%              Controls Group for Distributed and Uncertain Systems (CDUS)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

clear; clc; close all;

% Add all subfolders of the parent directory to the path
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')))

%% ====================== Data Setup ======================

orbit_type = 'NRHO_9_2';  % Specify orbit type (DRO, NRHO_L2_S, Halo_L1_N)
center     = 'bary';
p = BCR4BPOrbit(orbit_type, center);

Nrevs = 1;          % Number of revolutions to propagate
Nsteps = 10000;

order = 2;
scheme = 2;
%% ====================== Propagate Nrevs ======================

% Define propagation parameters
t0 = 0;                      % Initial time
tf = Nrevs * p.Tp;           % Final time = Nrevs full orbital periods
dt = p.Tp / Nsteps;          % Step size

% -------------------- ODE --------------------
opts = odeset('RelTol', 2e-13, 'AbsTol', 1e-13); %odeset('RelTol', 1e-16, 'AbsTol', 1e-30);
p.DS.integrator = @ode45;
tic
sol_ode = p.DS.propagate(p.nu0,[t0 tf],opts,@(t,x)p.DS.Hamiltons_EOM(t,x));
toc
sol_ode.t = sol_ode.x;
sol_ode.x = sol_ode.y;
sol_ode = rmfield(sol_ode,'y');
ODE_obj = Integrator(p,'RKF45', sol_ode);

% -------------------- SYMPLECTIC INTEGRATOR --------------------
% Create instance of the symplectic integrator class
SI_obj = SI(p, order, scheme);  

% Propagate using symplectic integrator
tic
SI_obj.propagate(t0, tf, dt);
toc
X_SI = SI_obj.sol.x; t_SI = SI_obj.sol.t;

% ------------------------ RUNGE-KUTTA ------------------------
% Create instance of the Runge-Kutta integrator class
RK_obj = RK(p, order);

% Propagate using RK integrator
tic
RK_obj.propagate(p.nu0, t0, tf, dt, ...
                                @(t, x) p.DS.Hamiltons_EOM(t, x));
toc
X_RK = RK_obj.sol.x; t_RK = RK_obj.sol.t;

%% ======================== POST-PROCESSING =========================

% ------------------------ Plot Orbits -----------------------------
ODE_obj.plot_traj(plot_2d_xy=true)                           
SI_obj.plot_traj(fig = gcf,plot_2d_xy=true)                           
RK_obj.plot_traj(fig = gcf,plot_2d_xy=true)

% ------------------- Plot Jacobi Constant Drift -------------------
% Plot the absolute Jacobi Constant difference from the initial value.
ODE_obj.plot_conserved_quantity(quantity='hamiltonian')                           
SI_obj.plot_conserved_quantity(fig = gcf,quantity='hamiltonian')
RK_obj.plot_conserved_quantity(fig = gcf,quantity='hamiltonian')