function u = util(param, c)
% CRRA utility u(c) = c^(1-gamma)/(1-gamma), gamma = param.nu (log if
% gamma = 1). Consumption is floored at param.minValue.
c = max(c, param.minValue);
if param.nu ~= 1
    u = c.^(1-param.nu) ./ (1-param.nu);
else
    u = log(c);
end
end