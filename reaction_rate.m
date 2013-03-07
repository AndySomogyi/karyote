function [w] = reaction_rate(r, y)
% Determines the rate of reaction for species i in compartment a. 
%
% r: reverse rate, n_comp x n_reactions
% y: affinity, n_comp x n_reactions
% returns w, size(n_comp, n_reactions)
    global R T
    w = r .* (exp(y / (R * T)) - 1);
end 

