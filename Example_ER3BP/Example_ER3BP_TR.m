%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Example Script for Symplectic Integrations
% 
% This script compares symplectic (SI) and Runge-Kutta (RK) integrators by
% propagating an orbit in the Circular Restricted Three-Body Problem
% (CR3BP) over several revolutions. It plots the 3D trajectory and the 
% Jacobi constant variation over time to evaluate energy preservation.
%
% Author: Anabel Soria-Carro 
% Date:   June 26, 2025
% Affiliation: The University of Texas at Austin
%              Controls Group for Distributed and Uncertain Systems (CDUS)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

clear; clc; close all;

% Add all subfolders of the parent directory to the path
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')))

%% ====================== Data Setup ======================

orbit_type = 'Halo_L2_3_1';  % Specify orbit type (DRO_1_1, DRO_3_2, DRO_5_2, Halo_L2_2_1, Halo_L2_3_1)
center     = 'p2';  % 'bary', 'p2'
Nrevs = 1;       % Number of revolutions to propagate

p = ER3BPOrbit(orbit_type, center, Nrevs);

order  = 4;         % Integrators order
scheme = 2;         % Störmer-Verlet scheme 1 or 2

reg_method = 'Russell'; % 'Sundman' or 'Russell'
%% ====================== Propagate Nrevs ======================

% -------------------- SYMPLECTIC INTEGRATOR --------------------
% Create instance of the symplectic integrator class
SI_obj = SI(p, order, scheme);

% Define propagation parameters
t0 = 0;                      % Initial time
tf = Nrevs * p.Tp;           % Final time = Nrevs full orbital periods
h = 0.04;%0.05;                 % Step size

TR_SI = TimeRegularized(SI_obj,reg_method);
TR_SI.propagate(t0, tf, h);

TR_SI.plot_traj_with_drift(font_size = 14,quantity='kamiltonian')

% ------------------------ RUNGE-KUTTA ------------------------
% Create instance of the Runge-Kutta integrator class
RK_obj = RK(p, order);

% Propagate using RK integrator
TR_RK = TimeRegularized(RK_obj,reg_method);
TR_RK.propagate(t0, tf, h);

TR_RK.plot_traj_with_drift(font_size = 14,quantity='kamiltonian')

%% -------------------- ODE --------------------
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12); 
p.DS.integrator = @ode45;
ODE_obj = Integrator(p,'RKF45');
tic
sol_ode = p.DS.propagate(p.nu0,[t0 tf],opts,@(t,x)p.DS.Hamiltons_EOM(t,x));
ODE_obj.time_solver = toc;
sol_ode.t = sol_ode.x;
sol_ode.x = sol_ode.y;
sol_ode = rmfield(sol_ode,'y');
sol_ode.Nsteps = round(length(sol_ode.t)/Nrevs);
sol_ode.coord = "hamiltonian";
ODE_obj.sol = sol_ode;

%% ======================== POST-PROCESSING =========================

% ------------------------ Plot Orbits -----------------------------
TR_SI.plot_traj(plot_2d_xy=true,show_steps=true)                           
TR_RK.plot_traj(fig = gcf,plot_2d_xy=true,show_steps=true)
ODE_obj.plot_traj(fig = gcf,plot_2d_xy=true,show_steps=true)                           

% ------------------- Plot Jacobi Constant Drift -------------------
% Plot the absolute Jacobi Constant difference from the initial value.
TR_SI.plot_conserved_quantity(quantity='kamiltonian',show_steps=true)
TR_RK.plot_conserved_quantity(fig = gcf,quantity='kamiltonian',show_steps=true)
ODE_obj.plot_conserved_quantity(fig = gcf,quantity='kamiltonian',show_steps=true)
