function Solution = DDCMSolver(DATA, qdata, Fdata)
%DDCMSolver Data-driven forward dynamics for Duffing oscillator (phase-space VI + ADM).
%
%   Solution = DDCMSolver(DATA, qdata, Fdata)
%
% Inputs
%   DATA  : struct with fields
%       m, c              : mass and damping
%       Dt, Tspan         : time step and [t0 tf]
%       Alpha             : symmetric Lagrangian parameter (0<=Alpha<=1)
%       Const             : reference stiffness Kref used in the metric
%       ForceFcn          : function handle Fext = ForceFcn(t,DATA)
%       Tol, TolC         : tolerances (equilibrium residual, cost change)
%       Iter              : max ADM iterations per time step
%       DetectJump        : iteration count used by a heuristic index-jump detector
%       p                 : metric exponent parameter forwarded to phiES
%       CI                : initial conditions [q0 qdot0]
%       Verbose (optional): true/false
%
%   qdata : cell array {eData_stage1, eData_stage2}
%   Fdata : cell array {fData_stage1, fData_stage2}
%
% Output (Solution struct)
%   Time      : time grid
%   Z         : [q; p] trajectory (2 x nSteps)
%   qmid,Fmid : continuous mid-point constitutive values at each step (2 x (nSteps-1))
%   qtil,Ftil : selected discrete constitutive pairs (2 x (nSteps-1))
%   Mu        : Lagrange multipliers from KKT solve (4 x nSteps)
%   Res       : equilibrium residual per step
%   Error     : infinity-norm change in KKT solution per ADM iteration (last iteration)
%   ErrorC    : infinity-norm change in cost per ADM iteration (last iteration)
%   Costq     : e-part cost
%   CostF     : f-part cost
%   Index     : cell array with index history for the selected data points
%   Jump      : per-step boolean indicating whether the jump heuristic activated
%   IndxJump  : indices of time steps with Jump==1
%
% Notes
%   - This function keeps the original numerical operations intact where possible,
%     but removes avoidable inv() calls and adds preallocation and documentation.
%   - The "jump" heuristic is preserved as in the original code.

arguments
    DATA (1,1) struct
    qdata (1,:) cell
    Fdata (1,:) cell
end

if ~isfield(DATA,'Verbose'); DATA.Verbose = false; end

% Build KKT blocks (constant for this problem)
[M, N, An, An1] = localKKTSystem(DATA);
KKT = [M, An.'; An, zeros(4,4)];

% Metric matrix CC = diag([Kref, Kref]) and its inverse (diagonal)
Kref = DATA.Const(1);
CC    = diag([Kref, Kref]);
invCC = diag([1/Kref, 1/Kref]);

% Time grid
h = DATA.Dt;
Time = DATA.Tspan(1):h:DATA.Tspan(2);
nSteps = numel(Time);
Alpha  = DATA.Alpha;

% Initial conditions in phase space
q0 = DATA.CI(1);
p0 = DATA.m * DATA.CI(2);
Z0 = [q0; p0];

% Initial discrete constitutive guess (consistent with the analytic spring law)
qtil_1a0 = q0;
qtil_a0  = q0;
Ftil_1a0 = DATA.KL*qtil_1a0 + DATA.KNL*qtil_1a0.^3;
Ftil_a0  = DATA.KL*qtil_a0  + DATA.KNL*qtil_a0.^3;
Y0       = [qtil_1a0; qtil_a0; Ftil_1a0; Ftil_a0];

% Preallocate outputs
Solution = struct();
Solution.Time  = Time;
Solution.Z     = zeros(2, nSteps);
Solution.Z(:,1)= Z0;
Solution.Mu    = zeros(4, nSteps);
Solution.qmid  = zeros(2, nSteps-1);
Solution.Fmid  = zeros(2, nSteps-1);
Solution.qtil  = zeros(2, nSteps-1);
Solution.Ftil  = zeros(2, nSteps-1);
Solution.Res   = zeros(1, nSteps-1);
Solution.Error = zeros(1, nSteps-1);
Solution.ErrorC= zeros(1, nSteps-1);
Solution.Costq = zeros(1, nSteps-1);
Solution.CostF = zeros(1, nSteps-1);
Solution.Iter  = zeros(1, nSteps-1);
Solution.Jump  = false(1, nSteps-1);
Solution.Index = cell(1, nSteps-1);

% Main time stepping
for i = 1:(nSteps-1)
    % Reset per-step quantities
    resid    = 1;
    errDX    = 1;
    errCost  = 1;
    Xprev    = zeros(10,1);
    costPrev = 0;
    IndexHist = [];
    jumpFlag = false;

    % Stage times (kept as in original implementation)
    t1a = (1-Alpha)*(i-1)*h + Alpha*i*h;
    ta  = Alpha*(i-1)*h     + (1-Alpha)*i*h;

    F1a = DATA.ForceFcn(t1a, DATA);
    Fa  = DATA.ForceFcn(ta,  DATA);

    k = 1;
    while ((resid >= DATA.Tol) || (errCost > DATA.TolC)) && (k <= DATA.Iter)
        % Continuous (equilibrium/compatibility) subproblem via KKT solve
        AUX = An1*Z0 + [ (h^2/2)*((1-Alpha)*F1a + Alpha*Fa);
                        (h/2)*(F1a + Fa);
                        0;
                        0];
        RHS = [N*Y0; AUX];
        X   = KKT \ RHS;

        qm = X(3:4,1);
        Fm = X(5:6,1);

        % Discrete (data) projection: ADM step
        if k > DATA.DetectJump
            % Original heuristic: if index history oscillates, freeze projection.
            diffIdx = abs(IndexHist(1,k-1) - IndexHist(1,k-3));
            if diffIdx == 0 || jumpFlag
                qtilde = Y0(1:2);
                Ftilde = Y0(3:4);
                jumpFlag = true;
            else
                [qtilde, Ftilde] = phiESOperator(qm, Fm, qdata, Fdata, DATA.p);
            end
        else
            [qtilde, Ftilde] = phiESOperator(qm, Fm, qdata, Fdata, DATA.p);
        end

        % Record selected indices (phiES returns exact data points)
        idx1 = find(qtilde(1) == qdata{1}, 1, 'first');
        idx2 = find(qtilde(2) == qdata{2}, 1, 'first');
        IndexHist = [IndexHist, [idx1; idx2]]; %#ok<AGROW>

        % Cost function (avoid inv(CC); CC is diagonal)
        Solution.Costq(i) = 0.5 * (qm - qtilde).' * CC    * (qm - qtilde);
        Solution.CostF(i) = 0.5 * (Fm - Ftilde).' * invCC * (Fm - Ftilde);
        cost = Solution.Costq(i) + Solution.CostF(i);

        % Residual of equilibrium constraint block
        resid = norm(An*X(1:6,1) - AUX, 1);

        % Convergence monitors
        errDX   = norm(X - Xprev, inf);
        errCost = norm(cost - costPrev, inf);

        % Update fixed-point variables
        Xprev    = X;
        costPrev = cost;
        Y0       = [qtilde; Ftilde];

        k = k + 1;
    end

    % Store per-step solution
    Solution.Z(:,i+1)   = X(1:2,1);
    Solution.qmid(:,i)  = X(3:4,1);
    Solution.Fmid(:,i)  = X(5:6,1);
    Solution.Mu(:,i+1)  = X(7:end,1);
    Solution.qtil(:,i)  = Y0(1:2);
    Solution.Ftil(:,i)  = Y0(3:4);
    Solution.Index{i}   = IndexHist;
    Solution.Jump(i)    = jumpFlag;
    Solution.Res(i)     = resid;
    Solution.Error(i)   = errDX;
    Solution.ErrorC(i)  = errCost;
    Solution.Iter(i)    = min(k-1, DATA.Iter);

    % Advance state
    Z0 = Solution.Z(:,i+1);

    if DATA.Verbose
        fprintf('Time step %d / %d | ADM iter = %d | Cost = %.7g | Res = %.3e\n', ...
            i, nSteps-1, Solution.Iter(i), costPrev, resid);
    end
end

Solution.IndxJump = find(Solution.Jump);

end

%% Local helpers
function [M, N, An, An1] = localKKTSystem(DATA)
% Constructs the constant blocks of the KKT system for the phase-space update.

h     = DATA.Dt;
m     = DATA.m;
c     = DATA.c;
Alpha = DATA.Alpha;
K     = DATA.Const(1);

An      = zeros(4,6);
An(1,1) = 1 + (h/(2*m))*c;
An(1,5) = (h^2/2)*((1-Alpha)/m);
An(1,6) = (h^2/2)*(Alpha/m);
An(2,1) = c;
An(2,2) = 1;
An(2,5) = h/2;
An(2,6) = h/2;
An(3,1) = -Alpha;
An(3,3) = 1;
An(4,1) = -(1-Alpha);
An(4,4) = 1;

An1      = zeros(4,2);
An1(1,1) = 1 + (h/(2*m))*c;
An1(1,2) = h/m;
An1(2,1) = c;
An1(2,2) = 1;
An1(3,1) = (1-Alpha);
An1(4,1) = Alpha;

M      = zeros(6,6);
M(3,3) = K;
M(4,4) = K;
M(5,5) = 1/K;
M(6,6) = 1/K;

N      = zeros(6,4);
N(3,1) = K;
N(4,2) = K;
N(5,3) = 1/K;
N(6,4) = 1/K;

end
