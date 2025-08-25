%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Example Script for Symplectic Integrations
% 
% This script compares symplectic (SI) and Runge-Kutta (RK) integrators by
% propagating an orbit in the Bicircular Restricted Four-Body Problem (BCR4BP)
% over several revolutions. 
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

orbit_type = 'NRHO_3_1';  % Specify orbit type (NRHO_3_1, NRHO_9_2)
center     = 'p2';
Nrevs = 1;                % Number of revolutions to propagate

p = BCR4BPOrbit(orbit_type, center,Nrevs);

order = 4;
scheme = 2;
%% ====================== Propagate Nrevs ======================

% Define propagation parameters
t0 = 0;                      % Initial time
tf = Nrevs * p.Tp;           % Final time = Nrevs full orbital periods
epsilon = 1;                 % Step size

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
ODE_obj.sol = sol_ode;

% -------------------- SYMPLECTIC INTEGRATOR --------------------
% Create instance of the symplectic integrator class
SI_obj = SI(p, order, scheme);  

% Propagate using symplectic integrator
params.alpha = 3/2;     % example value for time-regularization parameter

TR_SI = TimeRegularized(SI_obj,'sundman',params);
TR_SI.propagate(t0, tf, epsilon);

% ------------------------ RUNGE-KUTTA ------------------------
% Create instance of the Runge-Kutta integrator class
RK_obj = RK(p, order);

TR_RK = TimeRegularized(RK_obj,'sundman',params);
TR_RK.propagate(t0, tf, epsilon);

%TR_RK.plot_traj_with_drift(font_size = 14,quantity='hamiltonian')


%% ======================== POST-PROCESSING =========================

% ------------------------ Plot Orbits -----------------------------
ODE_obj.plot_traj(plot_2d_xy=false)                           
TR_SI.plot_traj(fig = gcf,plot_2d_xy=false)                           
TR_RK.plot_traj(fig = gcf,plot_2d_xy=false)

% ------------------- Plot Jacobi Constant Drift -------------------
% Plot the absolute Jacobi Constant difference from the initial value.
ODE_obj.plot_conserved_quantity(quantity='hamiltonian')                           
TR_SI.plot_conserved_quantity(fig = gcf,quantity='hamiltonian')
TR_RK.plot_conserved_quantity(fig = gcf,quantity='hamiltonian')