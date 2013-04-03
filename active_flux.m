function [ j_active ] = active_flux(c, s, k)
% Foo
%   c: concentration size: (single column vector of n_species * n_comp)
%   s: stochiometry matrix: ((n_species * n_comp), n_reactions)
%   k: reaction rate constants: n_reactions x 2

    n_reactions = length(k); 

    nu = zeros(n_reactions, 1);
    r = zeros(size(c));
    
    for i=1:n_reactions
        ind = find(s(:,i) < 0);
        r(:) = 0;
        r(ind) = abs(s(ind,i));
        
        nu(i) = k(i) * prod(c.^r);
    end
    
    j_active = s * nu;
end