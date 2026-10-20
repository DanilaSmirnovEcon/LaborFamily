function out = fiscal_welfare(agentPanel, param)
%FISCAL_WELFARE  Fiscal accounts and aggregate welfare of a simulated economy.
%   out = fiscal_welfare(agentPanel, param) returns per-capita monthly
%   averages over the simulated cross-sections (paper, Section 5):
%     out.T            net tax revenue: T_k(y) = y - D_k(y), credits netted
%     out.credits      tax credits paid: C_k(y) = D_k(y) - N(y)
%     out.B            unemployment benefits paid: b(a) of the unemployed
%     out.expenditure  fiscal expenditure: B + credits
%     out.balance      fiscal balance: T - B
%     out.W            aggregate welfare: population average of current
%                      (per-person) values
%   Taxable income y is the wage when employed and b(a) when unemployed;
%   home production is not taxed. Couples file jointly: T_k and C_k are
%   evaluated once on the pooled income and split equally between the
%   spouses, so each couple counts once.
%
%   Welfare changes in consumption-equivalent units, with CRRA gamma =
%   param.nu: CE = (W_policy / W_baseline)^(1/(1-gamma)) - 1.

emp  = agentPanel.employed == 1;
olf  = agentPanel.olf == 1;
unem = ~emp & ~olf;
mar  = agentPanel.married == 1;
kid  = agentPanel.has_kids == 1;

a = agentPanel.ability;
y = zeros(size(a));
y(emp)  = agentPanel.earnings(emp);
y(unem) = param.b(a(unem));
B = zeros(size(a));
B(unem) = param.b(a(unem));

% Spouse's taxable income (partner fields are NaN for singles)
pemp  = mar & (agentPanel.partner_employed == 1);
polf  = mar & (agentPanel.partner_olf == 1);
punem = mar & ~pemp & ~polf;
pa = agentPanel.partner_ability;
yp = zeros(size(a));
yp(pemp)  = agentPanel.partner_earnings(pemp);
yp(punem) = param.b(pa(punem));

yHH   = y + yp;                       % household taxable income
share = 1 - 0.5*mar;                  % per-person share of household taxes

N = param.tax(yHH, param);
D = N;
D(kid)  = param.onekid_tax(yHH(kid),  [], param);
D(~kid) = param.nokid_tax(yHH(~kid), [], param);

out.T           = mean(share(:) .* (yHH(:) - D(:)));
out.credits     = mean(share(:) .* (D(:) - N(:)));
out.B           = mean(B(:));
out.expenditure = out.B + out.credits;
out.balance     = out.T - out.B;
out.W           = mean(agentPanel.value(:), 'omitnan');
end
