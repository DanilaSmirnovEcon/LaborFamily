function tabs = policy_tables(results, param, outfile)
%POLICY_TABLES  Tables 7-10 of the paper from the policy counterfactuals.
%   tabs = policy_tables(results, param) takes the results struct built by
%   PolicyCounterfactuals.m (fields baseline, ub_plus15, credit_plus20,
%   each with moments_endog/moments_exog from compute_moments and
%   fw_endog/fw_exog from fiscal_welfare) and prints, relative to the
%   baseline:
%     Tables 7 and 9:  % change in average wage and unemployment rate by
%                      group (endogenous-marriage model)
%     Tables 8 and 10: % change in the marriage rate, aggregate welfare (CE)
%                      and fiscal expenditure, with and without endogenous
%                      marital choices
%   The tables are also written to outfile (default 'policy_tables.csv').

if nargin < 3 || isempty(outfile), outfile = 'policy_tables.csv'; end
gam = param.nu;

groups = {'never_married', 'Singles'; 'married', 'Married'; 'women', 'Women'; ...
          'men', 'Men'; 'haschild', 'Agents with a child'; 'nochild', 'Agents without a child'};
policies = {'ub_plus15',     'Tables 7-8: 15% higher unemployment benefits (b_min)'; ...
            'credit_plus20', 'Tables 9-10: 20% higher income caps for child-based tax credits'};

pc = @(new, old) 100*(new/old - 1);
ce = @(Wnew, Wold) 100*((Wnew/Wold)^(1/(1-gam)) - 1);

base = results.baseline;
rows = cell(0,4);   % {Table, Row, Column, Value}
for p = 1:size(policies,1)
    pol = results.(policies{p,1});
    fprintf('\n%s\n', policies{p,2});

    % Wages and unemployment rates by group (endogenous marriage)
    fprintf('%-20s', '');
    fprintf('%12s', 'Singles', 'Married', 'Women', 'Men', 'With child', 'No child'); fprintf('\n');
    wg = zeros(1,size(groups,1)); un = wg;
    for g = 1:size(groups,1)
        wg(g) = pc(pol.moments_endog.stats_wages.(groups{g,1}), base.moments_endog.stats_wages.(groups{g,1}));
        un(g) = pc(pol.moments_endog.(groups{g,1}).unemployed, base.moments_endog.(groups{g,1}).unemployed);
        rows(end+1,:) = {policies{p,1}, 'Wage (%)', groups{g,2}, wg(g)};               %#ok<AGROW>
        rows(end+1,:) = {policies{p,1}, 'Unemployment rate (%)', groups{g,2}, un(g)};  %#ok<AGROW>
    end
    fprintf('%-20s', 'Wage'); fprintf('%+11.2f%%', wg); fprintf('\n');
    fprintf('%-20s', 'Unemployment rate'); fprintf('%+11.2f%%', un); fprintf('\n\n');

    % Marriage rates, welfare and fiscal expenditure
    mr = pc(pol.moments_endog.married.all, base.moments_endog.married.all);
    wn = ce(pol.fw_endog.W, base.fw_endog.W);
    wx = ce(pol.fw_exog.W,  base.fw_exog.W);
    fn = pc(pol.fw_endog.expenditure, base.fw_endog.expenditure);
    fx = pc(pol.fw_exog.expenditure,  base.fw_exog.expenditure);
    fprintf('%-36s %15s %18s %20s\n', '', 'Marriage Rates', 'Aggregate Welfare', 'Fiscal Expenditure');
    fprintf('%-36s %+14.2f%% %+17.2f%% %+19.2f%%\n', 'With Endogenous Marital Choices', mr, wn, fn);
    fprintf('%-36s %15s %+17.2f%% %+19.2f%%\n', 'Without Endogenous Marital Choices', '-', wx, fx);
    rows(end+1,:) = {policies{p,1}, 'With endogenous marriage', 'Marriage rates (%)',        mr}; %#ok<AGROW>
    rows(end+1,:) = {policies{p,1}, 'With endogenous marriage', 'Aggregate welfare (CE, %)', wn}; %#ok<AGROW>
    rows(end+1,:) = {policies{p,1}, 'With endogenous marriage', 'Fiscal expenditure (%)',    fn}; %#ok<AGROW>
    rows(end+1,:) = {policies{p,1}, 'Without endogenous marriage', 'Aggregate welfare (CE, %)', wx}; %#ok<AGROW>
    rows(end+1,:) = {policies{p,1}, 'Without endogenous marriage', 'Fiscal expenditure (%)',    fx}; %#ok<AGROW>
end

fid = fopen(outfile, 'w');
fprintf(fid, 'Policy,Row,Column,Value\n');
for r = 1:size(rows,1)
    fprintf(fid, '"%s","%s","%s",%.10g\n', rows{r,1}, rows{r,2}, rows{r,3}, rows{r,4});
end
fclose(fid);
fprintf('\nSaved to %s\n', outfile);

tabs = [{'Policy','Row','Column','Value'}; rows];
end
