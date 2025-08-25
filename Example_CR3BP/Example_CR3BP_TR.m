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

orbit_type = 'NRHO_L2_S'; % Specify orbit type (DRO, NRHO_L2_S, Halo_L1_N)
center     = 'p2';  % 'bary', 'p2'
Nrevs = 1;          % Number of revolutions to propagate

p = CR3BPOrbit(orbit_type,center, Nrevs);

order  = 4;         % Integrators order
scheme = 2;         % Störmer-Verlet scheme 1 or 2

%% ====================== Propagate Nrevs ======================

% -------------------- SYMPLECTIC INTEGRATOR --------------------
% Create instance of the symplectic integrator class
SI_obj = SI(p, order, scheme);  

% Define propagation parameters
t0 = 0;                      % Initial time
tf = Nrevs * p.Tp;           % Final time = Nrevs full orbital periods
epsilon = 1;                 % Step size

% Propagate using symplectic integrator
params.alpha = 3/2;     % example value for time-regularization parameter

TR_SI = TimeRegularized(SI_obj,'sundman',params);
TR_SI.propagate(t0, tf, epsilon);

TR_SI.plot_traj_with_drift(font_size = 14,quantity='jacobi')

% ------------------------ RUNGE-KUTTA ------------------------
% Create instance of the Runge-Kutta integrator class
RK_obj = RK(p, order);

% Propagate using RK integrator
TR_RK = TimeRegularized(RK_obj,'sundman',params);
TR_RK.propagate(t0, tf, epsilon);

TR_RK.plot_traj_with_drift(font_size = 14,quantity='jacobi')

%% ======================== POST-PROCESSING =========================

% ------------------------ Plot Orbits -----------------------------
TR_SI.plot_traj(plot_2d_xy=true,show_steps=true)                           
TR_RK.plot_traj(fig = gcf,plot_2d_xy=true,show_steps=true)

% ------------------- Plot Jacobi Constant Drift -------------------
% Plot the absolute Jacobi Constant difference from the initial value.
TR_SI.plot_conserved_quantity(quantity='jacobi',show_steps=true)
TR_RK.plot_conserved_quantity(fig = gcf,quantity='jacobi',show_steps=true)
