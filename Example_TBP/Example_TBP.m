%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Example Script for Symplectic Integrations
% 
% This script compares symplectic (SI) and Runge-Kutta (RK) integrators by
% propagating an orbit in the Two-Body Problem (TBP) over several 
% revolutions. It plots the 3D trajectory and the energy variation over 
% time to evaluate energy preservation.
%
% Author: Anabel Soria-Carro 
% Date:   May 22, 2025
% Affiliation: The University of Texas at Austin
%              Controls Group for Distributed and Uncertain Systems (CDUS)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

clear; clc; close all;

% Add all subfolders of the parent directory to the path
addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'..')))

%% ====================== Data Setup ======================

orbit_type = 2;       % 1 or 2
main_body  = 'Earth';
p = TBPOrbit(main_body,orbit_type);

Nrevs  = 10;         % Number of revolutions to propagate
Nsteps = 1000;       % Number of time steps per revolution

order  = 8;         % Integrators order
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
[X_SI, t_SI] = SI_obj.propagate(t0, tf, dt);


% ------------------------ RUNGE-KUTTA ------------------------
% Create instance of the Runge-Kutta integrator class
RK_obj = RK(p, order);

% Propagate using RK integrator
[X_RK, t_RK] = RK_obj.propagate(p.nu0, t0, tf, dt, ...
                                @(t, x) p.DS.EOM(t, x));

%% ======================== POST-PROCESSING =========================

% ------------------------ Plot Orbits -----------------------------
SI_obj.plot_traj(plot_2d_xy=true)                           
RK_obj.plot_traj(fig = gcf,plot_2d_xy=true)

% ------------------- Plot Energy Drift -------------------
% Compute the energy over time for both integrators
% and plot the absolute difference from the initial value.

SI_obj.plot_conserved_quantity(quantity='energy')
RK_obj.plot_conserved_quantity(fig = gcf,quantity='energy')                               
