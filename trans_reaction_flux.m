function [flux] = trans_reaction_flux(c, s, k, z, v)
% trans compartment reaction rate via mass action.
% c: concentration size: (single column vector of n_species * n_comp)
%    organized by compartment, then species, 
%    C1S1, C1S2, C1S3, C2S1, C2S2,...
% s: stochiometry matrix: ((n_species * n_comp), n_reactions)
% k: reaction rate constants: n_reactions x 2
% z: valences, 1,n_species
% v: compartment voltages, n_comp, 1
%
% The flux from one compartment to the other, J, is given by:
% J = kC(1)
% where k is a constant for a given temperature and C is the concentration 
% in the compartment from which the flux is occurring.

    z = zeros(size(z));
    v = zeros(size(v));
    n_species = length(z);
    n_comp = length(v);
    flux = zeros(n_species, n_comp, n_comp);
    
    n_comp = length(v);
    n_reactions = size(k, 1); 
    nu = zeros(n_reactions, 1);
    
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
        rev = k(i,2) * prod(c.^rs(:,i));
        nu(i) = fwd - rev;
    end
    
    % matrix product of stochiometric matrix and rate vector, 
    % this is dc/dt for the intra-compartment part of the reactions. 
    dcdt = s * nu;
    
    s = reshape(s, [n_species, n_comp, n_reactions]);
    
    % heinously ugly, there's got to be a cleaner and faster
    % way of doing this, but its late...
    for k=1:n_reactions
        for i=1:n_comp
            outi = logical(s(:,i,k) < 0);
            % only look at other comps if we have outgoing flux.
            if any(outi)   
                for j=1:n_comp
                    if j ~= i
                        ini = logical(s(:,j,k) > 0);
                        fluxi = outi & ini;
                        flux(fluxi, i, j) = s(fluxi,j,k) * nu(k);
                    else
                    end
                end
            end
        end
    end
end 


