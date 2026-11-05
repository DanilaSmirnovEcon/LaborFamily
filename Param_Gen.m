function param = Param_Gen(policy)
%==========================================================================
% PARAM_GEN - Baseline parameters common to both genders
%==========================================================================
% Parameters of Morazzoni and Smirnov, "Labor and Family Dynamics in a
% Joint-Search Framework" (Table 2), in the paper's notation. Time is in
% months. Gender-specific parameters are set in Param_F.m and Param_M.m.
%
% Input (optional):
%   policy - function handle param = policy(param) applied to the scalar
%            parameters before the functions b(a), kappa(a), N(y), C_k(y)
%            are built (see PolicyCounterfactuals.m). Changing a scalar
%            after Param_Gen returns does NOT change those functions.
%
% Units: wages, benefits, child costs and home production are monthly
% flows in model consumption units (c.u.). The tax normalization
% ybar = param.averagewage = 61 c.u. corresponds to $45,000 a year
% ($25 an hour) in 2018 dollars: param.adjustW converts model wages to
% $/hour and param.adjustWScale converts annual dollars to c.u.
%==========================================================================

%==========================================================================
%% PREFERENCES AND DEMOGRAPHICS
%==========================================================================
param.nu   = 2.00;              % gamma: coefficient of relative risk aversion
param.rho  = 0.0034;            % rho: discount rate
param.zeta = 1/480;             % xi: retirement rate (40-year working life)

%==========================================================================
%% PRODUCTIVITY (paper eq. 1; common to both genders)
%==========================================================================
% employed:     da = phi_e  (abar  - a) dt + eps dB
% not employed: da = phi_ne (a_min - a) dt + eps dB
param.a.min = 0.00;             % a_min
param.a.max = 3.00;             % a_max
param.a.n   = 5;                % productivity grid points

param.phi_e    = 0.007;         % phi_e: productivity accumulation
param.phi_ne   = 2*param.phi_e; % phi_ne: productivity depletion (twice phi_e)
param.epsilon  = 0.065;         % eps: volatility of productivity shocks
param.abar     = 1.25;          % abar: long-run mean productivity

% Entry: a0 ~ logN(mu_a0, sigma_a0) truncated to [a_min, a_max]
param.mu_a0    = 0;             % mu_a0: normalized to 0
param.sigma_a0 = 0.365;         % sigma_a0: initial productivity dispersion

%==========================================================================
%% WAGE OFFERS
%==========================================================================
% Offers w ~ generalized Pareto (shape 0, scale sigma_w, location 0) on
% [w_min, w_max]. Realized earnings omega = theta_g w e^a, theta_m = 1.
param.w.min   = 10;             % w_min
param.w.max   = 30;             % w_max
param.w.n     = 9;              % wage grid points
param.sigma_w = 10;             % sigma_w
param.theta_w = 0.90;           % theta_w: gender wage-offer gap (women)

%==========================================================================
%% TAXES AND TRANSFERS
%==========================================================================
% N(y) = max{eta1 (y/ybar)^(1-eta2) ybar, 0.5 y}  (Heathcote et al., 2017)
param.lambda_tax  = 0.906;      % eta1: scale
param.tau_tax     = 0.181;      % eta2: progressivity
param.averagewage = 61;         % ybar: income normalization (c.u.)

param.adjustW      = param.averagewage/25;     % c.u. per $/hour
param.adjustWScale = param.averagewage/45000;  % c.u. per annual $

% Earned income tax credit C^M_k(y, l) (2018 EITC, paper eq. eitc_credit):
% phases in with household earnings l, phases out with household taxable
% income y (earnings plus unemployment benefits) above a threshold that
% depends on marital status M. Maximum credit Cbar_k = beta^Y1_k y^Y1_k.
param.pi1_n0 = 0.0765;          % beta^Y1_0: phase-in rate, no child
param.pi1_n1 = 0.34;            % beta^Y1_1: phase-in rate, child
param.pi2_n0 = 0.0765;          % beta^Y2_0: phase-out rate, no child
param.pi2_n1 = 0.1598;          % beta^Y2_1: phase-out rate, child
param.E1_n0  = 6750  * param.adjustWScale;   % y^Y1_0: end of phase-in, no child
param.E1_n1  = 10150 * param.adjustWScale;   % y^Y1_1: end of phase-in, child
param.E2s_n0 = 8490  * param.adjustWScale;   % y^Y2,s_0: start of phase-out, single, no child
param.E2s_n1 = 18660 * param.adjustWScale;   % y^Y2,s_1: start of phase-out, single, child
param.E2c_n0 = 14170 * param.adjustWScale;   % y^Y2,c_0: start of phase-out, couple, no child
param.E2c_n1 = 24350 * param.adjustWScale;   % y^Y2,c_1: start of phase-out, couple, child
% The credit reaches zero at y^Y2,M_k + Cbar_k/beta^Y2_k.

% Unemployment benefits b(a) = min{bmin + bscale e^a, bmax}
param.bmin   = 15.85;           % b_min (c.u.)
param.bscale = 0.85;            % b_scale
param.bmax   = 0.5*param.averagewage;   % b_max: 50% of the average wage

%==========================================================================
%% MARRIAGE AND FERTILITY
%==========================================================================
param.m   = 0.375;              % mu: marriage meeting rate
param.pcs = 0.0004;             % pi^s: fertility rate, singles
param.pcp = 0.005;              % pi^c: fertility rate, couples

% Child cost kappa(a) = kappa_min + kappa_scale e^a; couples pay one cost
% at the spouses' average productivity
param.cmin   = 7.0;             % kappa_min (c.u.)
param.cscale = 0.35;            % kappa_scale

% Forced marriage at a meeting (no value comparison); true only in the
% exogenous-marriage counterfactual (SimulatePanelCounter.m)
param.always_marry = false;

%==========================================================================
%% NUMERICAL SOLUTION PARAMETERS
%==========================================================================
% Trembling hand in the couples' switching decisions (paper Appendix B.2)
param.tremble = 1e-1;

% Value function iteration
param.updMult   = 0.15;         % update weight, value functions
param.updMult_s = 0.15;         % update weight, single-agent distributions
param.updMult_d = 0.5;          % update weight, couple distributions

% Convergence
param.maxiter = 100;            % maximum outer iterations
param.crit    = 1e-4;           % general convergence criterion
param.critC   = 1e-6;           % strict convergence criterion

% Grid refinement for single agents
param.Sdimmult = 1;
param.minValue = 1e-5;          % consumption floor in util.m

%==========================================================================
%% POLICY CHANGES (optional)
%==========================================================================
if nargin > 0 && ~isempty(policy)
    param = policy(param);
end

%==========================================================================
%% FUNCTIONS OF THE PARAMETERS (built after any policy change)
%==========================================================================
param.b = @(x) min(param.bmin + param.bscale.*exp(x), param.bmax);   % b(a)
param.c = @(x) param.cmin + param.cscale.*exp(x);                    % kappa(a)

% N(y): income after taxes, before credits
param.tax = @(x, param) max((param.lambda_tax .* ...
    ((x./param.averagewage).^(1-param.tau_tax)) .* param.averagewage), 0.5.*x);

% C^M_k(y, l): earned income tax credit; k = 1 with a child, M = 1 for a
% couple (filing jointly), 0 for a single
param.credit = @(y, l, k, M, param) eitc(y, l, k, M, param);

% D^M_k(y, l) = N(y) + C^M_k(y, l): disposable resources
param.D = @(y, l, k, M, param) param.tax(y, param) + param.credit(y, l, k, M, param);

%==========================================================================
%% WAGE OFFER DISTRIBUTION
%==========================================================================
% Probability of each wage grid node (men's offers; women's grid is the
% same distribution scaled by theta_w, see BuildSingleGrid/BuildDoubleGrid)
wvec = linspace(param.w.min, param.w.max, param.w.n);
dw   = (param.w.max - param.w.min)/(param.w.n - 1);
param.f = gppdf(wvec, 0, param.sigma_w, 0) .* dw;
param.f = param.f ./ sum(param.f);

% Same on the (optionally refined) single-agent wage grid
wvec_long = linspace(param.w.min, param.w.max, (param.w.n-1)*param.Sdimmult+1);
dw_long   = (param.w.max - param.w.min)/((param.w.n-1)*param.Sdimmult);
param.f_long = gppdf(wvec_long, 0, param.sigma_w, 0) .* dw_long;
param.f_long = param.f_long ./ sum(param.f_long);

end

%==========================================================================
%% EARNED INCOME TAX CREDIT C^M_k(y, l)
%==========================================================================

function C = eitc(y, l, k, M, param)
    % min{ beta1 l, Cbar, [Cbar - beta2 (y - y^Y2,M)^+]^+ }, with k and M
    % scalars; no credit without earnings
    if k
        b1 = param.pi1_n1;  b2 = param.pi2_n1;  y1 = param.E1_n1;
        if M, y2 = param.E2c_n1; else, y2 = param.E2s_n1; end
    else
        b1 = param.pi1_n0;  b2 = param.pi2_n0;  y1 = param.E1_n0;
        if M, y2 = param.E2c_n0; else, y2 = param.E2s_n0; end
    end
    Cbar = b1 * y1;
    C = min(min(b1 .* l, Cbar), max(Cbar - b2 .* max(y - y2, 0), 0));
end