function [gam, rho, info] = calibrate_exog_marriage(mode, targets, solutionSP_x, valueFunc_x, ...
                                                   gridSf, gridSm, gridP, param_x, simParam, opts)
%CALIBRATE_EXOG_MARRIAGE  Selection and sorting for the exogenous-marriage counterfactual.
%   [gam, rho, info] = calibrate_exog_marriage(mode, targets, solutionSP_x,
%       valueFunc_x, gridSf, gridSm, gridP, param_x, simParam)
%   chooses the selection (gam = param.exog_select) and sorting
%   (rho = param.exog_sort) parameters of SimulatePanelCounter so that the
%   economy without marital choice reproduces moments of the baseline with
%   endogenous marriage:
%     mode 'sorting':           gam = 0; rho matches targets.corr, the
%                               correlation between spouses' productivity
%                               (or targets.corr_field, e.g. 'corr_wages')
%     mode 'selection_sorting': gam matches targets.gap, the mean
%                               productivity of married minus singles
%                               (compute_moments: mean_ability.gap), and
%                               rho matches targets.corr
%   The marriage rate stays at param_x.exog_marriage_prob. Each evaluation
%   simulates the economy (SimulatePanelCounter with the same seed), so the
%   targets are matched on the same cross-sectional moments as the
%   baseline. Inputs as in PolicyCounterfactuals.m (values solved with
%   mu = 0). info records every evaluation.
%
%   opts (optional): tol_corr (default 0.01), tol_gap (0.01), maxit (8,
%   evaluations per one-dimensional search), rounds (3, alternating
%   searches in 'selection_sorting').

if nargin < 10, opts = struct(); end
if ~isfield(opts,'tol_corr'), opts.tol_corr = 0.01; end
if ~isfield(opts,'tol_gap'),  opts.tol_gap  = 0.01; end
if ~isfield(opts,'maxit'),    opts.maxit    = 8;    end
if ~isfield(opts,'rounds'),   opts.rounds   = 3;    end
if ~isfield(targets,'corr_field'), targets.corr_field = 'corr_ability'; end
hiBound = 0.99;

info.history = zeros(0,4);   % [gam, rho, corr, gap]

switch mode
    case 'sorting'
        gam = 0;
        rho = solve1(@(r) evalx(0, r, 'corr'), targets.corr, 0, hiBound, opts.tol_corr, opts.maxit);
    case 'selection_sorting'
        gam = 0;  rho = 0;
        for r = 1:opts.rounds
            gam = solve1(@(g) evalx(g, rho, 'gap'), targets.gap, 0, hiBound, opts.tol_gap, opts.maxit);
            rho = solve1(@(x) evalx(gam, x, 'corr'), targets.corr, 0, hiBound, opts.tol_corr, opts.maxit);
            [c, g] = evalx(gam, rho, 'corr');
            if abs(c - targets.corr) < opts.tol_corr && abs(g - targets.gap) < opts.tol_gap
                break
            end
        end
    otherwise
        error('calibrate_exog_marriage: unknown mode %s', mode);
end

[info.corr, info.gap] = evalx(gam, rho, 'corr');
info.targets = targets;
fprintf('calibrate_exog_marriage (%s): gamma = %.3f, rho = %.3f; %s = %.3f (target %.3f), gap = %.3f (target %s)\n', ...
    mode, gam, rho, targets.corr_field, info.corr, targets.corr, info.gap, num2str(getfield_or(targets,'gap',NaN)));

    function [c, g] = evalx(gm, rh, which)
        % Simulated moments at (gm, rh); repeated points are taken from history
        hit = find(abs(info.history(:,1) - gm) < 1e-12 & abs(info.history(:,2) - rh) < 1e-12, 1);
        if ~isempty(hit)
            c = info.history(hit,3);  g = info.history(hit,4);
        else
            p = param_x;  p.exog_select = gm;  p.exog_sort = rh;
            panel = SimulatePanelCounter(solutionSP_x, valueFunc_x, gridSf, gridSm, gridP, p, simParam);
            m = compute_moments(panel, simParam, p);
            c = m.(targets.corr_field);  g = m.mean_ability.gap;
            info.history(end+1,:) = [gm, rh, c, g];
        end
        if strcmp(which, 'gap'), [c, g] = deal(g, c); end   % return the searched moment first
    end
end

function x = solve1(f, target, lo, hi, tol, maxit)
% Root of f(x) = target on [lo, hi] by the Illinois (regula falsi) method;
% f is increasing in x. If the target is out of reach, the closer bound is
% returned with a warning.
    a = lo;  fa = f(a) - target;
    if abs(fa) < tol, x = a; return; end
    b = hi;  fb = f(b) - target;
    if abs(fb) < tol, x = b; return; end
    if sign(fa) == sign(fb)
        if abs(fa) < abs(fb), x = a; else, x = b; end
        warning('calibrate_exog_marriage:outOfReach', ...
            'Target %.3f not reached on [%.2f, %.2f]; using %.2f.', target, lo, hi, x);
        return
    end
    side = 0;
    for it = 1:maxit
        x  = (a*fb - b*fa) / (fb - fa);
        fx = f(x) - target;
        if abs(fx) < tol, return; end
        if sign(fx) == sign(fb)
            b = x;  fb = fx;
            if side == -1, fa = fa/2; end
            side = -1;
        else
            a = x;  fa = fx;
            if side == 1, fb = fb/2; end
            side = 1;
        end
    end
end

function v = getfield_or(s, f, d)
    if isfield(s, f), v = s.(f); else, v = d; end
end
