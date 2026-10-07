function T = calibration_table(m, param, outfile)
%CALIBRATION_TABLE  Model moments next to the paper's data targets.
%   T = calibration_table(moments, param) prints every moment reported in
%   the paper (Tables 3, 4, 5 and the couples' employment shares) with the
%   model value, the data value and the gap. moments is the output of
%   compute_moments. T is a cell array with a header row:
%   {Section, Moment, Model, Data, Gap, Previous}.
%
%   The table is also written to outfile (default 'calibration_table.csv').
%   If that file already exists from a previous run, its model values are
%   shown in a 'Previous' column, so the effect of a parameter change can
%   be read off directly.
%
%   Units follow the paper: rates and shares in %, wages in 2018 dollars.

if nargin < 3 || isempty(outfile), outfile = 'calibration_table.csv'; end
if isfield(param,'adjustW') && ~isempty(param.adjustW), adjW = param.adjustW; else, adjW = 1; end

rows = cell(0,4);   % {Section, Moment, Model, Target}
    function add(section, name, model, target)
        rows(end+1,:) = {section, name, model, target};
    end
pct = @(path) 100*getf(m, path);

% ---- Table 3: targeted moments ----
s = 'Targeted';
add(s, 'EU rate, men (%)',                          pct('transitions.male.EU'),     1.5);
add(s, 'EU rate, women (%)',                        pct('transitions.female.EU'),   1.2);
add(s, 'UE rate, men (%)',                          pct('transitions.male.UE'),     26);
add(s, 'UE rate, women (%)',                        pct('transitions.female.UE'),   23);
add(s, 'OE rate, men (%)',                          pct('transitions.male.OE'),     5.2);
add(s, 'OE rate, women (%)',                        pct('transitions.female.OE'),   4.0);
add(s, 'Share 1-yr wage change < 20% (%)',          pct('wagechangeyear'),          65);
add(s, 'Wage change per extra month unemployed (%)', pct('wagechangeunemployment_per_month_logpct'), -1);
add(s, 'UB as % of income (%)',                     pct('unemp_benefits.all_share_of_median_income'), 35);
add(s, 'Delta UB, bottom vs top earners (p.p.)', ...
    100*(getf(m,'unemp_benefits.lowinc_mean_share') - getf(m,'unemp_benefits.highinc_mean_share')), 15);
add(s, 'Gender wage gap (%)', ...
    100*(1 - getf(m,'stats_wages.women')/getf(m,'stats_wages.men')), 20);
add(s, 'Wage progression (entry/average)',          getf(m,'wageprogression'),       0.55);
add(s, 'Wage dispersion (var log wage)',            getf(m,'wagedisp'),              0.70);
add(s, 'Initial wage dispersion (var log wage)',    getf(m,'wageentryvar'),          0.35);
add(s, 'Dispersion in job offers (CV)',             getf(m,'stats_wages.cvOffers'),  0.40);
add(s, 'Share of married (%)',                      pct('married.all'),              52);
add(s, 'Child cost as % of income (%)',             pct('child_cost.mean_share'),    15);
add(s, 'Delta child cost, bottom vs top (p.p.)', ...
    100*(getf(m,'child_cost.low.mean_share') - getf(m,'child_cost.high.mean_share')), 15);
add(s, 'OLF single women, no kids (%)',             pct('olf_cells.single_women_nokid'),   12.18);
add(s, 'OLF single women, kids (%)',                pct('olf_cells.single_women_kid'),     19.99);
add(s, 'OLF single men, no kids (%)',               pct('olf_cells.single_men_nokid'),     10.71);
add(s, 'OLF single men, kids (%)',                  pct('olf_cells.single_men_kid'),        8.01);
add(s, 'OLF married women, no kids (%)',            pct('olf_cells.married_women_nokid'),  14.51);
add(s, 'OLF married women, kids (%)',               pct('olf_cells.married_women_kid'),    25.82);
add(s, 'OLF married men, no kids (%)',              pct('olf_cells.married_men_nokid'),     3.59);
add(s, 'OLF married men, kids (%)',                 pct('olf_cells.married_men_kid'),       3.08);

% ---- Table 4: untargeted labor-market and marriage moments ----
s = 'Untargeted';
add(s, 'Avg hourly wage, men ($)',        getf(m,'stats_wages.men')/adjW,           26.5);
add(s, 'Avg hourly wage, women ($)',      getf(m,'stats_wages.women')/adjW,         21.6);
add(s, 'Avg hourly wage, married ($)',    getf(m,'stats_wages.married')/adjW,       26.0);
add(s, 'Avg hourly wage, singles ($)',    getf(m,'stats_wages.never_married')/adjW, 21.2);
add(s, 'Avg hourly wage, with kids ($)',  getf(m,'stats_wages.haschild')/adjW,      25.0);
add(s, 'Avg hourly wage, no kids ($)',    getf(m,'stats_wages.nochild')/adjW,       21.7);
add(s, 'Unemployment rate, men (%)',      pct('men.unemployed'),           5.09);
add(s, 'Unemployment rate, women (%)',    pct('women.unemployed'),         3.72);
add(s, 'Unemployment rate, married (%)',  pct('married.unemployed'),       3.44);
add(s, 'Unemployment rate, singles (%)',  pct('never_married.unemployed'), 6.51);
add(s, 'Unemployment rate, with kids (%)', pct('haschild.unemployed'),     3.90);
add(s, 'Unemployment rate, no kids (%)',  pct('nochild.unemployed'),       5.28);
add(s, 'Wage correlation between spouses', getf(m,'corr_wages'),           0.35);
add(s, 'Age at first marriage, men',      getf(m,'age_at_marriage.men_years'),   28.3);
add(s, 'Age at first marriage, women',    getf(m,'age_at_marriage.women_years'), 25.9);

% ---- Table 5: wage marital premium and unemployment marital gap ----
s = 'Table 5';
add(s, 'Wage marital premium (%)', ...
    100*(getf(m,'stats_wages.married')/getf(m,'stats_wages.never_married') - 1), 20.83);
add(s, 'Unemployment marital gap (p.p.)', ...
    100*(getf(m,'never_married.unemployed') - getf(m,'married.unemployed')), 3.07);

% ---- Shares of couples by employment status ----
s = 'Couples';
add(s, 'Couples (E,E) (%)',       pct('couplesemploymentshares.EE'), 55.0);
add(s, 'Couples (E,O)+(O,E) (%)', pct('couplesemploymentshares.EO'), 32.4);
add(s, 'Couples (E,U)+(U,E) (%)', pct('couplesemploymentshares.EU'),  4.2);
add(s, 'Couples (U,U) (%)',       pct('couplesemploymentshares.UU'),  0.4);
add(s, 'Couples (U,O)+(O,U) (%)', pct('couplesemploymentshares.UO'),  0.8);
add(s, 'Couples (O,O) (%)',       pct('couplesemploymentshares.OO'),  8.2);

model  = cell2mat(rows(:,3));
target = cell2mat(rows(:,4));
gap    = model - target;

% Previous run, matched by moment name
prev = NaN(size(model));
if exist(outfile, 'file')
    [pNames, pModel] = read_previous(outfile);
    [found, loc] = ismember(rows(:,2), pNames);
    prev(found) = pModel(loc(found));
end

% Print
fprintf('\nMODEL VS. DATA (paper Tables 3-5)\n');
fprintf('%-46s %10s %10s %10s %10s\n', 'Moment', 'Model', 'Data', 'Gap', 'Previous');
section = '';
for r = 1:size(rows,1)
    if ~strcmp(rows{r,1}, section)
        section = rows{r,1};
        fprintf('--- %s ---\n', section);
    end
    fprintf('%-46s %10.3f %10.3f %10.3f %10.3f\n', rows{r,2}, model(r), target(r), gap(r), prev(r));
end

% Save (moment names contain commas, so text fields are quoted)
fid = fopen(outfile, 'w');
fprintf(fid, 'Section,Moment,Model,Data,Gap\n');
for r = 1:size(rows,1)
    fprintf(fid, '"%s","%s",%.10g,%.10g,%.10g\n', rows{r,1}, rows{r,2}, model(r), target(r), gap(r));
end
fclose(fid);

T = [{'Section','Moment','Model','Data','Gap','Previous'}; ...
     rows(:,1:2), num2cell([model, target, gap, prev])];
fprintf('Saved to %s\n', outfile);
end

function [names, model] = read_previous(outfile)
% Moment names and model values from a CSV written by this function.
names = {}; model = [];
fid = fopen(outfile, 'r');
if fid < 0, return; end
fgetl(fid);                                   % header
line = fgetl(fid);
while ischar(line)
    tok = regexp(line, '^"[^"]*","([^"]*)",([^,]*),', 'tokens', 'once');
    if numel(tok) == 2
        names{end+1,1} = tok{1};              %#ok<AGROW>
        model(end+1,1) = str2double(tok{2});  %#ok<AGROW>
    end
    line = fgetl(fid);
end
fclose(fid);
end

function v = getf(s, path)
% Value of a nested field (e.g. 'transitions.male.EU'), or NaN if missing.
v = NaN;
try
    parts = strsplit(path, '.');
    for p = 1:numel(parts), s = s.(parts{p}); end
    if isnumeric(s) && isscalar(s), v = double(s); end
catch
end
end
