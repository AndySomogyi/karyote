function [w] = intra_reaction_rate(c, s, k)
% Intra compartment reaction rate via chemical potential. 
%   c: concentration size: (single column vector of n_species * n_comp)
%   s: stochiometry matrix: ((n_species * n_comp), n_reactions)
%   k: reaction rate constants: n_reactions x 2

    n_reactions = size(k, 1); 
    rev = zeros(n_reactions, 1);
    
    % chemical potential without reference bit
    mu = log(c);
    
    % reference chem potential
    g0 = log(k(:,1) ./ k(:,2));
    
    a = (g0 - s' * mu);
    
    % reverse part of stochiometry matrix
    rs = s;
    rs(rs < 0) = 0;
    
    fs = s;
    fs(fs > 0) = 0;
    fs = abs(fs);
    
    for i = 1:n_reactions
        rev(i) = k(i,2) * prod(c.^rs(:,i));
    end
    
    nu = rev .* (exp(a) - 1);
    
    w = s * nu;
    
end 

