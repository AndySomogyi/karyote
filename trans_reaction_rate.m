function [rate] = trans_reaction_rate(c, s, k, n_comp)
% Trans compartment reaction rate via mass action
%   c: concentration size: (single column vector of n_species * n_comp)
%   s: stochiometry matrix: ((n_species * n_comp), n_reactions)
%   k: reaction rate constants: n_reactions x 2

    n_reactions = size(k, 1); 
    nu = zeros(n_reactions, 1);
    n_species = size(s, 1) / n_comp;
    flux = zeros(n_species,n_comp,n_comp);
    
   
    rs = s;
    rs(rs < 0) = 0;
    
    fs = s;
    fs(fs > 0) = 0;
    fs = abs(fs);
    
    for i = 1:n_reactions
        fwd = k(i,1) * prod(c.^fs(:,i));
        rev = k(i,2) * prod(c.^rs(:,i));
        nu(i) = fwd - rev;
    end
    
    rate = s * nu;
    
    for i=1:n_comp
        for j=i:n_comp
            fwd = s((i:i+n_species),:);
            
            
        end
    end
end 


