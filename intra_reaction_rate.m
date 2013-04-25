function [dcdt] = intra_reaction_rate(c, s, k, z, v)
% Intra compartment reaction rate via chemical potential. 
% c: concentration size: (single column vector of n_species * n_comp)
%    organized by compartment, then species, 
%    C1S1, C1S2, C1S3, C2S1, C2S2,...
% s: stochiometry matrix: ((n_species * n_comp), n_reactions)
% k: reaction rate constants: n_reactions x 2
% z: valences, 1,n_species
% v: compartment voltages, n_comp, 1

    global F R T
    
    n_comp = length(v);
    
    n_reactions = size(k, 1); 
    nu = zeros(n_reactions, 1);
    n_species = size(s, 1) / n_comp;
    
    % potential for all species in all compartments
    v = kron(v, ones(n_species, 1));
    
    % valence for all species in all compartments, 
    % make it a column vector (TODO, this should start as a column vec).
    z = kron(ones(1, n_comp), z)';
    
    % electro part of chemical potential, 
    mu_e = ((F/(R*T)) * z .* v)';
   
    % reverse stoichiometric constants
    rs = s;
    rs(rs < 0) = 0;
    
    % foward stochiometric constants
    fs = s;
    fs(fs > 0) = 0;
    fs = abs(fs);
    
    % could be one line each, but breaking up makes debugging easier. 
    for i = 1:n_reactions
        fwd = k(i,1) * prod(c.^fs(:,i));
        fwd_e = mu_e * fs(:,i) / (R*T);
        rev = k(i,2) * prod(c.^rs(:,i));
        rev_e = mu_e * rs(:,i) / (R*T);
        nu(i) = fwd * exp(fwd_e) - rev * exp(rev_e);
    end
    
    % matrix product of stochiometric matrix and rate vector, 
    % this is dc/dt for the intra-compartment part of the reactions. 
    dcdt = s * nu;
end 

