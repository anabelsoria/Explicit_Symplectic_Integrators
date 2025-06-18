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

Nrevs = 10;          % Number of revolutions to propagate
Nsteps = 1000;       % Number of time steps per revolution

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
X_SI = SI_obj.sol.x; t_SI = SI_obj.sol.t;

SI_obj.plot_traj_dCfluct(title="SI" + num2str(order))

% ------------------------ RUNGE-KUTTA ------------------------
% Create instance of the Runge-Kutta integrator class
RK_obj = RK(p, order);

% Propagate using RK integrator
tic
RK_obj.propagate(p.nu0, t0, tf, dt, ...
                                @(t, x) p.DS.Hamiltons_EOM(t, x));
toc
X_RK = RK_obj.sol.x; t_RK = RK_obj.sol.t;

RK_obj.plot_traj_dCfluct(title="RK" + num2str(order))

%% ======================== POST-PROCESSING =========================

% ------------------------ Plot Orbits -----------------------------
SI_obj.plot_traj(plot_2d_xy=true)                           
RK_obj.plot_traj(fig = gcf,plot_2d_xy=true)

% ------------------- Plot Jacobi Constant Drift -------------------
% Plot the absolute Jacobi Constant difference from the initial value.
SI_obj.plot_dC
RK_obj.plot_dC(fig = gcf)
