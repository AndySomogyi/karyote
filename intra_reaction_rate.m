function [dcdt] = intra_reaction_rate(c, s, k, z, v)
% Intra compartment reaction rate via chemical potential. 
% c: concentration size: (single column vector of n_species * n_comp)
%    organized by compartment, then species, 
%    C1S1, C1S2, C1S3, C2S1, C2S2,...
% s: stochiometry matrix: ((n_species * n_comp), n_reactions)
% k: reaction rate constants: n_reactions x 2
% z: valences, m species
% v: compartment voltages, n comp

    global F R T

    n_reactions = size(k, 1); 
    rev = zeros(n_reactions, 1);
    
    % chemical potential without reference bit
    mu = log(c);
    
    % reference chem potential
    g0 = log(k(:,1) ./ k(:,2));
    
    % electrochemical part of chem potential, the Kronecker product
    % is used to unfold the potential and valences to the same layout
    % as the concentrations. 
    mu_e = reshape((F/(R*T))*kron(v, z), size(mu));
    
    mu = mu + mu_e;
    
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
    
    % the net reaction rate. 
    nu = rev .* (exp(a) - 1);
    
    % matrix product of stochiometric matrix and rate vector, 
    % this is dc/dt for the intra-compartment part of the reactions. 
    dcdt = s * nu;
end 

