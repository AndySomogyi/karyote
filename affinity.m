function [ a ] = affinity(g_delta, s, mu)
% The total affinity for the reactions. 
% adds up the saved reference chem pot, and the stoich weighted
% chemical potentials for each species in the reaction. 
%
% g_delta: sum of reference chem potential, (n_reaction, 1)
% s: stochiometry matrix: ((n_species * n_comp), n_reactions)
% mu: chemical poteial for each species, with out the reference chem pot. 
%     (n_species, 1)
%    

    
    a = s' * mu + g_delta;

end