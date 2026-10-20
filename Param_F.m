function [param, gridSf] = Param_F(param)
%==========================================================================
% PARAM_F - Parameters specific to female agents
%==========================================================================
% Sets the gender-specific parameters of the paper (Table 2) and builds
% the female single-agent grids.
%
% Inputs:
%   param  - Base parameter structure from Param_Gen()
%
% Outputs:
%   param  - Parameter structure with param.pf added
%   gridSf - State space grids for female single agents
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

pf.mu       = param.phi_e;      % employed: da = phi_e (abar - a) dt + eps dB
pf.mu_u     = param.phi_ne;     % not employed: da = phi_ne (a_min - a) dt + eps dB
pf.a_sigma  = param.epsilon;    % volatility of productivity shocks
pf.amean    = param.abar;       % long-run mean productivity when employed
pf.aergmean = param.mu_a0;      % entry: a0 ~ logN(mu_a0, sigma_a0), truncated
pf.aergsd   = param.sigma_a0;
pf.pcs      = param.pcs;        % fertility rate as a single (pi^s)

%==========================================================================
%% LABOR MARKET PARAMETERS
%==========================================================================

% Wage offers: theta_g = theta_w for women (paper Section 2.2)
pf.wscale = param.theta_w;

% Job separation (delta_g)
pf.delta = 0.015;

% Job arrival rates
pf.lambda_u  = 0.45;          % unemployed (lambda_g)
pf.lambda_oE = 0.105;         % out of the labor force (nu_g), direct O -> W
pf.lambda_o  = 50;            % switching intensity into/out of the labor force
                                % (numerical device, paper Appendix B.2)

%==========================================================================
%% HOME PRODUCTION (out of the labor force; paper Table 2, c.u.)
%==========================================================================
% Not taxable. Constant by marital status and parenthood.

pf.hs0 = 30.26;                  % single, no child  (h^s_{g,k=0})
pf.hs1 = 29.61;                  % single, child     (h^s_{g,k=1})
pf.hc0 = 27.11;                  % married, no child (h^c_{g,k=0})
pf.hc1 = 26.46;                  % married, child    (h^c_{g,k=1})

%==========================================================================
%% BUILD FEMALE SINGLE AGENT GRIDS
%==========================================================================

[gridSf, pf] = BuildSingleGrid(param, pf);

%==========================================================================
%% ENTRY DISTRIBUTION
%==========================================================================

% Productivity of newborn female agents (paper: a0 ~ log-normal truncated
% to [a_min, a_max]). Each grid node gets the probability of its own cell
% (midpoints between nodes, half-cells at the edges), renormalized for the
% truncation. SimulatePanel draws a0 from the same truncated distribution.
edgesA   = [gridSf.a.min; (gridSf.a.a(1:end-1) + gridSf.a.a(2:end))/2; gridSf.a.max];
pf.a.erg = diff(logncdf(edgesA, pf.aergmean, pf.aergsd));
pf.a.erg = pf.a.erg ./ sum(pf.a.erg);

%==========================================================================
%% STORE IN PARAMETER STRUCTURE
%==========================================================================

param.pf = pf;

end