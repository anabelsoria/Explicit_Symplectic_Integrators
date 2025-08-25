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

purpleMatlab = [0.4940 0.1840 0.5560];
yellowMatlab = [0.9290 0.6940 0.1250];
greenMatlab  = [0.4660 0.6740 0.1880];
blueMatlab = [0 0.4470 0.7410];
orangeMatlab = [0.8500 0.3250 0.0980];
%% ====================== Data Setup ======================

orbit_type = 'DPO'; % Specify orbit type (DRO, NRHO_L2_S, Halo_L1_N,Lyapunov_L2)
center     = 'bary';
Nrevs = 6;          % Number of revolutions to propagate

p = CR3BPOrbit(orbit_type, center, Nrevs);

Nsteps = 1000;      % Number of time steps per revolution

order  = 4;         % Integrators order
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
SI_obj.propagate(t0, tf, dt);
toc
X_SI = SI_obj.sol.x; t_SI = SI_obj.sol.t; SI_obj.sol.Nsteps = Nsteps;

SI_obj.plot_traj_with_drift(title="SI" + num2str(order),font_size = 14,quantity='jacobi')

% ------------------------ RUNGE-KUTTA ------------------------
% Create instance of the Runge-Kutta integrator class
RK_obj = RK(p, order);

% Propagate using RK integrator
tic
RK_obj.propagate(p.nu0, t0, tf, dt, ...
                                @(t, x) p.DS.Hamiltons_EOM(t, x));
toc
X_RK = RK_obj.sol.x; t_RK = RK_obj.sol.t; RK_obj.sol.Nsteps = Nsteps;

RK_obj.plot_traj_with_drift(title="RK" + num2str(order),font_size = 14,quantity='jacobi')

% -------------------- ODE --------------------
opts = odeset('RelTol', 2e-13, 'AbsTol', 1e-13); %odeset('RelTol', 1e-16, 'AbsTol', 1e-30);
p.DS.integrator = @ode45;
ODE_obj = Integrator(p,'RKF45');
tic
sol_ode = p.DS.propagate(p.nu0,[t0 tf],opts,@(t,x)p.DS.Hamiltons_EOM(t,x));
ODE_obj.time_solver = toc;
sol_ode.t = sol_ode.x;
sol_ode.x = sol_ode.y;
sol_ode = rmfield(sol_ode,'y');
sol_ode.nsteps = round(length(sol_ode.t)/Nrevs);
sol_ode.coord = 'hamiltonian';
sol_ode.Nsteps = round(length(sol_ode.x)/p.Nrevs);
ODE_obj.sol = sol_ode;

%% ======================== POST-PROCESSING =========================

% ------------------------ Plot Orbits -----------------------------
ODE_obj.plot_traj(plot_2d_xy=true,font_size = 14,color='b')       
SI_obj.plot_traj(fig = gcf,plot_2d_xy=true,font_size = 14,color=greenMatlab)                           
RK_obj.plot_traj(fig = gcf,plot_2d_xy=true,font_size = 14,color=orangeMatlab)

% ------------------- Plot Jacobi Constant Drift -------------------
% Plot the absolute Jacobi Constant difference from the initial value.
ODE_obj.plot_conserved_quantity(quantity='jacobi',font_size = 14,color='b',show_steps = true)                           
SI_obj.plot_conserved_quantity(fig = gcf,quantity='jacobi',font_size = 14,color=greenMatlab,show_steps = true)
RK_obj.plot_conserved_quantity(fig = gcf,quantity='jacobi',font_size = 14,color=orangeMatlab,show_steps = true)

%% ------------------- Plot Error with ODE -------------------
[~,ode_x] = p.DS.propagate(p.nu0,SI_obj.sol.t,opts,@(t,x)p.DS.Hamiltons_EOM(t,x));

SI_state_error = vecnorm(SI_obj.sol.x-ode_x');
SI_obj.plot_state_error(SI_state_error,font_size = 14,color=greenMatlab)

RK_state_error = vecnorm(RK_obj.sol.x-ode_x');
RK_obj.plot_state_error(RK_state_error,fig=gcf,font_size = 14,color=orangeMatlab)

