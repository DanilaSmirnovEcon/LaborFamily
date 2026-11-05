function tab = table5_decomposition(r, outfile)
%TABLE5_DECOMPOSITION  Table 5 with intermediate counterfactuals.
%   tab = table5_decomposition(r) takes results.baseline from
%   PolicyCounterfactuals.m and compares, for the marital wage premium and
%   the unemployment marital gap:
%     (1) without endogenous marriage, random matching  (moments_exog)
%     (2) without endogenous marriage, with sorting      (moments_exog_sort)
%     (3) without endogenous marriage, with selection
%         and sorting                                    (moments_exog_selsort)
%     (4) with endogenous marriage                       (moments_endog)
%   Columns (2) and (3) reproduce the baseline's correlation between
%   spouses' productivity, and (3) also its productivity gap between
%   married and single individuals (calibrate_exog_marriage.m). Moving from
%   (1) to (2) measures the role of sorting, from (2) to (3) of selection
%   into marriage, and from (3) to (4) of marital choice itself (timing of
%   marriage, the option to stay single, and the expectation of marriage).
%   The table is also written to outfile (default 'table5_decomposition.csv').

if nargin < 2 || isempty(outfile), outfile = 'table5_decomposition.csv'; end

cols  = {'moments_exog', 'moments_exog_sort', 'moments_exog_selsort', 'moments_endog'};
names = {'Random matching', 'Sorting', 'Selection+sorting', 'Endogenous'};
have  = isfield(r, cols);

rowsNames = {'Marital wage premium (%)', 'Unemployment marital gap (p.p.)', ...
             'Corr. spouses'' productivity', 'Corr. spouses'' wages', ...
             'Productivity gap, married - single', 'Share married (%)'};
data = [20.83, 3.07, NaN, 0.35, NaN, 52];

vals = NaN(numel(rowsNames), numel(cols));
for c = 1:numel(cols)
    if ~have(c), continue; end
    m = r.(cols{c});
    vals(:,c) = [100*(m.stats_wages.married/m.stats_wages.never_married - 1);
                 100*(m.never_married.unemployed - m.married.unemployed);
                 m.corr_ability;
                 m.corr_wages;
                 m.mean_ability.gap;
                 100*m.married.all];
end

fprintf('\nTABLE 5: LABOR MARKET OUTCOMES WITH AND WITHOUT ENDOGENOUS MARRIAGE\n');
fprintf('%-36s %9s', '', 'Data'); fprintf(' %18s', names{:}); fprintf('\n');
for i = 1:numel(rowsNames)
    fprintf('%-36s %9.2f', rowsNames{i}, data(i)); fprintf(' %18.2f', vals(i,:)); fprintf('\n');
end
if isfield(r, 'exog_calibration')
    k = r.exog_calibration;
    fprintf('Sorting: rho = %.3f; selection and sorting: gamma = %.3f, rho = %.3f\n', ...
        k.sort.rho, k.selsort.gam, k.selsort.rho);
end

fid = fopen(outfile, 'w');
fprintf(fid, 'Row,Data,%s\n', strjoin(names, ','));
for i = 1:numel(rowsNames)
    fprintf(fid, '"%s",%.10g', rowsNames{i}, data(i));
    fprintf(fid, ',%.10g', vals(i,:)); fprintf(fid, '\n');
end
fclose(fid);
fprintf('Saved to %s\n', outfile);

tab = [{'Row','Data'}, names; rowsNames(:), num2cell([data(:), vals])];
end
