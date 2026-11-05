%==========================================================================
% POLICYCOUNTERFACTUALS.M - Section 5 policy experiments (Tables 7-10)
%==========================================================================
% Policy counterfactuals of Morazzoni & Smirnov, "Labor and Family
% Dynamics in a Joint-Search Framework", Section 5:
%   1. ub_plus15      - 15% higher minimum unemployment benefit b_min
%   2. credit_plus20  - 20% higher income caps {y^Y1_1, y^Y2_1} for the
%                       child-based tax credit (households with a child,
%                       single and married)
% Other parameters stay at the baseline (Param_Gen, Param_F, Param_M).
%
% For each policy and the baseline the model is solved and simulated:
%   - with endogenous marriage: SolveEndogProblem + SimulatePanel
%   - without endogenous marriage (Table 5 counterfactual): values are
%     re-solved with a zero meeting rate (param.m = 0), so nobody expects
%     to meet, and SimulatePanelCounter marries 52% of newborns at birth:
%       exog         - to a random partner (random matching)
%       exog_sort    - with sorting: spouses' productivity correlation as
%                      in the baseline with endogenous marriage
%       exog_selsort - with selection and sorting: also the baseline's
%                      productivity gap between married and singles
%     The sorting and selection parameters are calibrated once, on the
%     baseline (calibrate_exog_marriage.m), and held fixed under the
%     policies, like the marriage rate.
%
% Outputs: Table 5 with the intermediate counterfactuals
% (table5_decomposition.m, table5_decomposition.csv) and Tables 7-10
% (policy_tables.m, policy_tables.csv), plus all moments, fiscal accounts
% and welfare in PolicyCounterfactuals_results.mat.
%
% Run manually (not called from Main.m). Set run_policies = false for the
% baseline only (Table 5). The calibration of the intermediate
% counterfactuals adds about 10-30 simulations of the baseline.
%==========================================================================

clear all; close all;

run_policies = true;    % false: baseline only (Table 5)
run_sorting  = true;    % intermediate counterfactuals (sorting, selection)
calibTargets.corr_field = 'corr_ability';   % or 'corr_wages'

% Policy changes are applied inside Param_Gen, before b(a), kappa(a) and
% the tax functions are built, so that the functions reflect them.
policies = struct( ...
    'name',  {'baseline',  'ub_plus15',        'credit_plus20'}, ...
    'apply', {[],          @apply_ub_plus15,   @apply_credit_plus20} ...
);
if ~run_policies, policies = policies(1); end

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

    % ---- 5. Store ----
    results.(name).param         = param;
    results.(name).moments_endog = moments_endog;
    results.(name).moments_exog  = moments_exog;
    results.(name).fw_endog      = fw_endog;
    results.(name).fw_exog       = fw_exog;

    % ---- 4c. Without endogenous marriage, with sorting / selection ----
    if run_sorting
        if strcmp(name, 'baseline')
            % Calibrate to the baseline with endogenous marriage
            calibTargets.corr = moments_endog.(calibTargets.corr_field);
            calibTargets.gap  = moments_endog.mean_ability.gap;
            [~, rS, infoS] = calibrate_exog_marriage('sorting', calibTargets, ...
                solutionSP_x, valueFunc_x, gridSf, gridSm, gridP, param_x, simParam);
            [gSS, rSS, infoSS] = calibrate_exog_marriage('selection_sorting', calibTargets, ...
                solutionSP_x, valueFunc_x, gridSf, gridSm, gridP, param_x, simParam);
            exogCal.sort    = struct('gam', 0,   'rho', rS,  'info', infoS);
            exogCal.selsort = struct('gam', gSS, 'rho', rSS, 'info', infoSS);
            results.(name).exog_calibration = exogCal;
        end
        variants = {'sort', 'selsort'};
        for v = 1:numel(variants)
            cv = exogCal.(variants{v});
            pv = param_x;  pv.exog_select = cv.gam;  pv.exog_sort = cv.rho;
            panel = SimulatePanelCounter(solutionSP_x, valueFunc_x, gridSf, gridSm, gridP, pv, simParam);
            results.(name).(['moments_exog_' variants{v}]) = compute_moments(panel, simParam, pv);
            results.(name).(['fw_exog_' variants{v}])      = fiscal_welfare(panel, pv);
        end
    end

    save(sprintf('solutionSP_%s.mat', name), 'solutionSP', 'param', 'gridSf', 'gridSm', 'gridP', '-v7.3');
end

save('PolicyCounterfactuals_results.mat', 'results', '-v7.3');

% ---- 6. Table 5 (with/without endogenous marriage) and Tables 7-10 ----
table5_decomposition(results.baseline);
if run_policies
    policy_tables(results, results.baseline.param);
end

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
    % households with a dependent child only (k = 1), regardless of marital
    % status: the end of the phase-in (y^Y1_1, which also raises the
    % maximum credit beta^Y1_1 y^Y1_1) and the start of the phase-out for
    % singles and couples (y^Y2,s_1 and y^Y2,c_1).
    param.E1_n1  = param.E1_n1  * 1.20;
    param.E2s_n1 = param.E2s_n1 * 1.20;
    param.E2c_n1 = param.E2c_n1 * 1.20;
end
