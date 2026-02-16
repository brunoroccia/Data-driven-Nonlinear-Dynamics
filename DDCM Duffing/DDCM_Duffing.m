%% DDCM_Duffing.m
% Data-Driven Computational Mechanics (DDCM) for the forced Duffing oscillator
% using an Alternating Direction Method (ADM) projection and a variational
% (phase-space) time integrator.
%
% This script is an executable example that:
%   1) defines the Duffing oscillator and time-integration parameters,
%   2) builds a 1D constitutive dataset (e,f),
%   3) calls DDCMSolver to compute the DDCM trajectory,
%   4) produces a constitutive-space diagnostic plot.
%
% Repository structure (expected on MATLAB path):
%   - DDCM_Duffing.m      (this file)
%   - DDCMSolver.m        (DDCM forward-dynamics solver)
%   - MyForce.m           (external forcing definition)
%   - phiESOperator.m     (ADM data projection wrapper)
%   - phiES.m             (single spring data projection)
%
% (c) Open-source release for journal publication.

clear; clc; close all;

%% Problem definition (Duffing oscillator)
%   m*qdd + c*qd + KL*q + KNL*q^3 = Fext(t)
DATA = struct();
DATA.m     = 2.0;      % Mass
DATA.c     = 0.5;      % Damping coefficient
DATA.KL    = 1.0;      % Linear stiffness
DATA.KNL   = 0.5;      % Cubic stiffness

% External forcing
DATA.F     = 1.0;      % Force amplitude
DATA.W     = 0.5;      % Forcing frequency
DATA.Force = 'MyForce';
DATA.Flag  = 0;        % 0: cosine
                      % 1: sine with smooth ramp on
                      % 2: sigmoid ramp (step-like)
                      % 3: smoothstep C^2 ramp

% Integrator / algorithm parameters
DATA.Alpha      = 0.5;     % Symmetric Lagrangian parameter (alpha=1/2 => time-symmetric)

% Time discretization
DATA.Dt    = 0.025;
DATA.Tspan = [0 200];

% ADM / fixed-point tolerances
DATA.Tol        = 1e-12;   % Equilibrium residual tolerance
DATA.TolC       = 1e-6;    % Cost change tolerance
DATA.Iter       = 50;      % Max ADM iterations per time step
DATA.DetectJump = 4;       % Heuristic: detect index "jumps" after this many iterations

% Metric parameter in phiES (see phiES.m)
DATA.p = 1;

% Initial conditions: [q0, qdot0]
DATA.CI = [1.0, 0.0];

% Optional runtime flags
DATA.Verbose = true;       % Print per-step iteration summary

%% Constitutive dataset generation (1D spring law)
% Dataset represents the nonlinear spring: f = KL*e + KNL*e^3
nd = 101;
qmin = -1.5;
qmax = +1.5;

eq = linspace(qmin, qmax, nd);
fq = DATA.KL*eq + DATA.KNL*eq.^3;

% For this Duffing example, the algorithm stores two identical datasets,
% associated with the (1-alpha) and alpha stages.
qdata = {eq, eq};
Fdata = {fq, fq};

% Build a reference linearized stiffness to form the diagonal metric matrix.
% (The solver uses this to weight the e- and f- distances.)
idxNonZero = (eq ~= 0);
Kref = mean(fq(idxNonZero) ./ eq(idxNonZero));
DATA.Const = Kref;

%% Solve
DATA.ForceFcn = str2func(DATA.Force);
Solution = DDCMSolver(DATA, qdata, Fdata);

%% Plot: constitutive space (dataset + selected discrete pairs)
figure('Name','DDCM Duffing: Constitutive space');
plot(eq, fq, 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'b');
hold on;
plot(Solution.qtil(1,:), Solution.Ftil(1,:), 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'k');
plot(Solution.qtil(2,:), Solution.Ftil(2,:), 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'r');
plot(Solution.qmid(1,:), Solution.Fmid(1,:), 'x', 'MarkerSize', 6, 'MarkerEdgeColor', 'c');

xlabel('Spring elongation e');
ylabel('Spring force f');
legend('Dataset', 'Selected pair for 1-\alpha', 'Selected pair for \alpha', 'Continuous (mid) for 1-\alpha', 'Location','Best');
grid on;
