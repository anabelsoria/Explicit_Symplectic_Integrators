%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Example Script for Symplectic Integrations
% 
% This script compares symplectic (SI) and Runge-Kutta (RK) integrators by
% propagating an orbit in the Elliptic Restricted Three-Body Problem
% (ER3BP) over several revolutions.
%
% Author: Anabel Soria-Carro 
% Date:   June 19, 2025
% Affiliation: The University of Texas at Austin
%              Controls Group for Distributed and Uncertain Systems (CDUS)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

clear; clc; close all;

% Add all subfolders of the parent directory to the path
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')))

%% ====================== Data Setup ======================

orbit_type = 'DRO';  % Specify orbit type (DRO)
center     = 'bary';
p = ER3BPOrbit(orbit_type, center);

Nrevs = 5E3;          % Number of revolutions to propagate
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
K_SI = SI_obj.prob.DS.kamiltonian(X_SI);

SI_obj.plot_traj_with_drift(title="SI" + num2str(order),font_size = 14,quantity='kamiltonian')

% ------------------------ RUNGE-KUTTA ------------------------
% Create instance of the Runge-Kutta integrator class
RK_obj = RK(p, order);

% Propagate using RK integrator
tic
RK_obj.propagate(p.nu0, t0, tf, dt, ...
                                @(t, x) p.DS.Hamiltons_EOM(t, x));
toc
X_RK = RK_obj.sol.x; t_RK = RK_obj.sol.t;
K_RK = SI_obj.prob.DS.kamiltonian(X_RK);

RK_obj.plot_traj_with_drift(title="RK" + num2str(order),font_size = 14,quantity='kamiltonian')

%% ======================== POST-PROCESSING =========================

% ------------------------ Plot Orbits -----------------------------
SI_obj.plot_traj(plot_2d_xy=true, font_size = 14)                           
RK_obj.plot_traj(fig = gcf,plot_2d_xy=true, font_size = 14)

% ------------------- Plot Kamiltonian Drift -------------------
% Plot the absolute Kamiltonian difference from the initial value.
SI_obj.plot_conserved_quantity(quantity='kamiltonian')
RK_obj.plot_conserved_quantity(fig = gcf,quantity='kamiltonian')
