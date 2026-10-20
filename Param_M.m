function [param, gridSm] = Param_M(param)
%==========================================================================
% PARAM_M - Parameters specific to male agents
%==========================================================================
% Sets the gender-specific parameters of the paper (Table 2) and builds
% the male single-agent grids.
%
% Inputs:
%   param  - Base parameter structure from Param_Gen()
%
% Outputs:
%   param  - Parameter structure with param.pm added
%   gridSm - State space grids for male single agents
%
% Gender-specific (paper Table 2): wage-offer wedge theta_g, job separation
% delta_g, job arrival rates lambda_g (unemployed) and nu_g (out of the
% labor force), home production h^s_{g,k} and h^c_{g,k}. Productivity,
% entry and fertility parameters are common to both genders (Param_Gen).
%
% Reference: Morazzoni, M. and Smirnov, D., "Labor and Family Dynamics in
% a Joint-Search Framework"
%==========================================================================

%==========================================================================
%% PRODUCTIVITY AND ENTRY (common to both genders, set in Param_Gen)
%==========================================================================

pm.mu       = param.phi_e;      % employed: da = phi_e (abar - a) dt + eps dB
pm.mu_u     = param.phi_ne;     % not employed: da = phi_ne (a_min - a) dt + eps dB
pm.a_sigma  = param.epsilon;    % volatility of productivity shocks
pm.amean    = param.abar;       % long-run mean productivity when employed
pm.aergmean = param.mu_a0;      % entry: a0 ~ logN(mu_a0, sigma_a0), truncated
pm.aergsd   = param.sigma_a0;
pm.pcs      = param.pcs;        % fertility rate as a single (pi^s)

%==========================================================================
%% LABOR MARKET PARAMETERS
%==========================================================================

% Wage offers: theta_g = 1 for men (paper Section 2.2)
pm.wscale = 1;

% Job separation (delta_g)
pm.delta = 0.012;

% Job arrival rates
pm.lambda_u  = 0.55;          % unemployed (lambda_g)
pm.lambda_oE = 0.095;         % out of the labor force (nu_g), direct O -> W
pm.lambda_o  = 50;            % switching intensity into/out of the labor force
                                % (numerical device, paper Appendix B.2)

%==========================================================================
%% HOME PRODUCTION (out of the labor force; paper Table 2, c.u.)
%==========================================================================
% Not taxable. Constant by marital status and parenthood.

pm.hs0 = 23.20;                  % single, no child  (h^s_{g,k=0})
pm.hs1 = 22.50;                  % single, child     (h^s_{g,k=1})
pm.hc0 = 20.05;                  % married, no child (h^c_{g,k=0})
pm.hc1 = 19.35;                  % married, child    (h^c_{g,k=1})

%==========================================================================
%% BUILD MALE SINGLE AGENT GRIDS
%==========================================================================

[gridSm, pm] = BuildSingleGrid(param, pm);

%==========================================================================
%% ENTRY DISTRIBUTION
%==========================================================================

% Productivity of newborn male agents (paper: a0 ~ log-normal truncated
% to [a_min, a_max]). Each grid node gets the probability of its own cell
% (midpoints between nodes, half-cells at the edges), renormalized for the
% truncation. SimulatePanel draws a0 from the same truncated distribution.
edgesA   = [gridSm.a.min; (gridSm.a.a(1:end-1) + gridSm.a.a(2:end))/2; gridSm.a.max];
pm.a.erg = diff(logncdf(edgesA, pm.aergmean, pm.aergsd));
pm.a.erg = pm.a.erg ./ sum(pm.a.erg);

%==========================================================================
%% STORE IN PARAMETER STRUCTURE
%==========================================================================

param.pm = pm;

end