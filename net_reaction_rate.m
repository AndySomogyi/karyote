function [w] = net_reaction_rate(k, a)
% Determines the rate of reaction for species i in compartment a. 
%
% y: affinity, n_comp x n_reactions
% returns w, size(n_comp, n_reactions)
    global R T
    w = k(:,2) .* (exp(a / (R * T)) - 1);
end 

