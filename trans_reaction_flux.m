function [flux] = trans_reaction_flux(c, s, k, ~, ~, o, ~)
% trans compartment reaction rate via mass action.
% c: concentration size: (single column vector of n_species * n_comp)
%    organized by compartment, then species, 
%    C1S1, C1S2, C1S3, C2S1, C2S2,...
% s: stochiometry matrix: ((n_species * n_comp), n_reactions)
% k: reaction rate constants: n_reactions x 2
% o: volume of each compartment, size: n_comp,1
    
    n_comp = size(o,1);
    n_species = size(s,1) / n_comp;

    % z: valences, 1,n_species
    z = zeros(1,n_species);
    % v: compartment voltages, n_comp, 1
    v = zeros(n_comp,1);

    % dc dt: change in concentration over time, mol/(volume*seconds), 
    % however the k constants are defined to be different units to 
    % actually yield a flux, flux: mol/(area*seconds)
    dcdt = intra_reaction_rate(c, s, k, z, v)
end 


