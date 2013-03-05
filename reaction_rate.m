

function [w] = reaction_rate(r, y)
% r: reverse rate, k species x n compartments
% y: affinity, k species x n compartments
    global R T
    w = r .* (exp(y / (R * T)) - 1);
end 

