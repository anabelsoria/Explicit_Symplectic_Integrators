%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Example Script for the SI Precision Ladder
%
% Propagates the same CR3BP orbit with the symplectic integrator's four
% arithmetic options (double, ICS, CS, dd) and compares their energy drift.
%
% Author: Anabel Soria-Carro
% Date:   September 3, 2026
% Affiliation: The University of Texas at Austin
%              Controls Group for Distributed and Uncertain Systems (CDUS)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

clear; clc; close all;

% Add all subfolders of the parent directory to the path
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')))

c = plot_colors();

%% ====================== Data Setup ======================

orbit_type = 'DRO'; % Specify orbit type (DRO, NRHO_L2_S, Halo_L1_N,Lyapunov_L2)
center     = 'p2';  % 'bary', 'p2'
Nrevs = 10;         % Number of revolutions to propagate. CS and dd carry
                     % real per-step cost, so this stays modest here; push
                     % it toward the 5000-10000 range of Example_CR3BP_TR.m
                     % to see the retention gap between precisions widen.

p = CR3BPOrbit(orbit_type, center, Nrevs);

Nsteps = 10000;       % Number of time steps per revolution

order  = 4;         % Integrators order
scheme = 2;         % Only Stormer-Verlet B carries the ICS/CS/dd kernels

t0 = 0;
tf = Nrevs * p.Tp;
dt = p.Tp / Nsteps;

%% ====================== Propagate Nrevs ======================

precisions = {'double','ics','cs','dd'};
colors     = {c.blue, c.orange, c.green, c.purple};

SI_objs = cell(size(precisions));
for k = 1:numel(precisions)
    SI_objs{k} = SI(p, order, scheme, precisions{k});
    SI_objs{k}.name = SI_objs{k}.name + " (" + precisions{k} + ")";
    tic
    SI_objs{k}.propagate(t0, tf, dt);
    toc
    SI_objs{k}.sol.Nsteps = Nsteps;
end

%% ======================== POST-PROCESSING =========================

% ------- Plot Jacobi Constant Fluctuation Averaged over Each Rev -------
fig = [];
for k = 1:numel(precisions)
    if isempty(fig)
        SI_objs{k}.plot_conserved_quantity(quantity='jacobi',font_size = 14, ...
            color=colors{k},avg_per_rev=true);
        fig = gcf;
    else
        SI_objs{k}.plot_conserved_quantity(fig=fig,quantity='jacobi', ...
            font_size = 14,color=colors{k},avg_per_rev=true);
    end
end

% ------------------------ Plot Orbits -----------------------------
fig = [];
for k = 1:numel(precisions)
    if isempty(fig)
        SI_objs{k}.plot_traj(plot_2d_xy=true,font_size = 14,color=colors{k});
        fig = gcf;
    else
        SI_objs{k}.plot_traj(fig=fig,plot_2d_xy=true,font_size = 14,color=colors{k});
    end
end

% -------------------- Plot Wall-Clock Time per Precision --------------------
times = cellfun(@(o) o.time_solver, SI_objs);

figure; hold on; grid on;
for k = 1:numel(precisions)
    bar(k, times(k), 'FaceColor', colors{k});
end
set(gca,'XTick',1:numel(precisions),'XTickLabel',upper(precisions));
xlabel('Precision', 'Interpreter', 'latex');
ylabel('Wall-clock time [s]', 'Interpreter', 'latex');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14);
