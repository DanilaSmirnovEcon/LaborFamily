%==========================================================================
% POLICYCOUNTERFACTUALS.M - Section 5 policy experiments (Tables 7-10)
%==========================================================================
% Policy counterfactuals of Morazzoni & Smirnov, "Labor and Family
% Dynamics in a Joint-Search Framework", Section 5:
%   1. ub_plus15      - 15% higher minimum unemployment benefit b_min
%   2. credit_plus20  - 20% higher income caps {y^Y1_1, y^Y2_1} for the
%                       child-based tax credit (households with a child)
% Other parameters stay at the baseline (Param_Gen, Param_F, Param_M).
%
% For each policy and the baseline the model is solved and simulated
% twice:
%   - endogenous marriage: SolveEndogProblem + SimulatePanel
%   - without endogenous marriage (Table 5 counterfactual): values are
%     re-solved with a zero meeting rate (param.m = 0), so nobody expects
%     to meet, and SimulatePanelCounter marries 52% of newborns to a
%     random partner at birth
%
% Outputs: Tables 7-10 printed and saved to policy_tables.csv
% (policy_tables.m), plus all moments, fiscal accounts and welfare in
% PolicyCounterfactuals_results.mat.
%
% Run manually (not called from Main.m); it takes about 6x as long as Main.
%==========================================================================

clear all; close all;

% Policy changes are applied inside Param_Gen, before b(a), kappa(a) and
% the tax functions are built, so that the functions reflect them.
policies = struct( ...
    'name',  {'baseline',  'ub_plus15',        'credit_plus20'}, ...
    'apply', {[],          @apply_ub_plus15,   @apply_credit_plus20} ...
);

simParam = struct( ...
  'N',                 1000, ...
  'seed',              12345, ...
  'nWorkers',          [], ...
  'coarse_periods',    69, ...
  'coarse_step_months',12, ...
  'burn_periods',      200, ...
  'burn_step_months',  0.25, ...
  'sim_periods',       480, ...
  'sim_step_months',   0.25 ...
);

results = struct();

for k = 1:numel(policies)
    name = policies(k).name;
    fprintf('\n==========================================================================\n');
    fprintf('Policy: %s\n', name);
    fprintf('==========================================================================\n');

    % ---- 1. Parameters and grids ----
    param = Param_Gen(policies(k).apply);
    [param, gridSf] = Param_F(param);
    [param, gridSm] = Param_M(param);
    [param, gridP]  = BuildDoubleGrid(param, true);
    param.NSall = param.pf.NSall + param.pm.NSall;

    % ---- 2. Single-agent problems ----
    param.updMult = 0.90;
    param.critC   = 1e-6;
    param.maxiter = 5000;
    solutionSf = SolveSingleProblem(param, param.pf, gridSf);
    solutionSm = SolveSingleProblem(param, param.pm, gridSm);

    % ---- 3. Endogenous marriage equilibrium ----
    param.critC     = 1e-6;
    param.updMult   = 0.75;
    param.updMult_s = 0.5;
    param.updMult_d = 0.75;
    param.maxiter   = 200;
    solutionSP = SolveEndogProblem(solutionSf, solutionSm, false, true, ...
                                   param, gridSf, gridSm, gridP);
    valueFunc  = approximateValueFunctions(solutionSP, gridSf, gridSm, gridP, param);

    % ---- 4a. Simulate with endogenous marriage ----
    agentPanel_endog = SimulatePanel(solutionSP, valueFunc, gridSf, gridSm, gridP, param, simParam);
    moments_endog = compute_moments(agentPanel_endog, simParam, param);
    fw_endog      = fiscal_welfare(agentPanel_endog, param);

    % ---- 4b. Without endogenous marriage: values with no meetings ----
    param_x   = param;
    param_x.m = 0;
    solutionSP_x = SolveEndogProblem(solutionSf, solutionSm, false, true, ...
                                     param_x, gridSf, gridSm, gridP);
    valueFunc_x  = approximateValueFunctions(solutionSP_x, gridSf, gridSm, gridP, param_x);
    agentPanel_exog = SimulatePanelCounter(solutionSP_x, valueFunc_x, gridSf, gridSm, gridP, param_x, simParam);
    moments_exog = compute_moments(agentPanel_exog, simParam, param_x);
    fw_exog      = fiscal_welfare(agentPanel_exog, param_x);

    % ---- 5. Store and save ----
    results.(name).param         = param;
    results.(name).moments_endog = moments_endog;
    results.(name).moments_exog  = moments_exog;
    results.(name).fw_endog      = fw_endog;
    results.(name).fw_exog       = fw_exog;

    save(sprintf('solutionSP_%s.mat', name), 'solutionSP', 'param', 'gridSf', 'gridSm', 'gridP', '-v7.3');
end

save('PolicyCounterfactuals_results.mat', 'results', '-v7.3');

% ---- 6. Table 5 (with/without endogenous marriage) and Tables 7-10 ----
mb = results.baseline.moments_endog;  mx = results.baseline.moments_exog;
fprintf('\nTable 5: marital wage premium / unemployment marital gap\n');
fprintf('  without endogenous marriage: %+6.2f%% / %+5.2f p.p.\n', ...
    100*(mx.stats_wages.married/mx.stats_wages.never_married - 1), ...
    100*(mx.married.unemployed - mx.never_married.unemployed));
fprintf('  with endogenous marriage:    %+6.2f%% / %+5.2f p.p.\n', ...
    100*(mb.stats_wages.married/mb.stats_wages.never_married - 1), ...
    100*(mb.married.unemployed - mb.never_married.unemployed));

policy_tables(results, results.baseline.param);

%==========================================================================
%% POLICY PARAMETER CHANGES (applied inside Param_Gen)
%==========================================================================

function param = apply_ub_plus15(param)
    % Table 7-8: b(a) = min{b_min + b_scale e^a, b_max} with a 15% higher
    % minimum b_min. Other benefit parameters unchanged.
    param.bmin = param.bmin * 1.15;
end

function param = apply_credit_plus20(param)
    % Table 9-10: 20% higher income caps for the child-based tax credit,
    % households with a dependent child only (k = 1): the end of the
    % phase-in (y^Y1_1) and the start of the phase-out (y^Y2_1).
    param.E1_n1 = param.E1_n1 * 1.20;
    param.E2_n1 = param.E2_n1 * 1.20;
end
