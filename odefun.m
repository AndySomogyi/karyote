function [ f ] = odefun(cap, a, l, h, z, o, s, k)
% Create the odefunc that is given to an integrator. 
% 
% Args: 
%    cap: compartment capacitance, size(n_comp,n_comp), 
%        must have zero diagonal. 
%    a: area between compartments (membrane), size(n_comp, n_comp), 
%        must have zero diagonal. 
%    l: thickness of membrane separating a and a' compartments, 
%        size(n_comp, n_comp), must have zero diagonal. 
%    h: membrane permeability for species i between compartments, 
%        size(n_species, n_comp, n_comp).
%    z: valence of species i, size(n_species)
%    o: volume of compartments, size(n_comp)
%    s: stochiometry matrix: (n_species, n_comp, n_reactions)
%    k: forward and backward rates, n_reaction x 2.


    f = @fun;
    
    % all sizes are determined by the stoichiometry matrix, everything
    % else should agree.
    n_species = size(s, 1);
    n_comp = size(s, 2);
    n_reactions = size(s, 3);
    
    s = reshape(s, [n_species*n_comp, n_reactions]);
    
    % total membrane flus of species i from a to a'
    j = zeros(n_species, n_comp, n_comp);
    
    % charge-neutral flux of species from a' to a, 
    jn = zeros(size(j));
    
    dcdt_trans = zeros(n_species, n_comp);
    
    % h: membrane permeability. size: n_species x n_comp x n_comp.
    assert(size(h, 1) == n_species, 'permeability must be [n_species,n_comp,n_comp]');
    assert(size(h, 2) == n_comp,    'permeability must be [n_species,n_comp,n_comp]');
    assert(size(h, 3) == n_comp,    'permeability must be [n_species,n_comp,n_comp]');
    
    
    function [state] = fun(t,state)
        % the entire state, [concentration; potential] is wrapped
        % up in the state vector s. s should be 
        % (n_comp*n_species) + n_comp long 
        
        % concentration
        % c: size(n_comp, n_comp)
        %
        % voltage
        % v: size(1,n_comp)
        
        fprintf('time: %d\n', t);
        
        % concentration vector, organized by compartment, then species. 
        c = state(1:(n_comp*n_species));
        % voltage vector, n_comp x 1
        v = state((n_comp*n_species)+1:end);
        
        e = equilibrium_factor(c,z,v,l);
        j(:) = membrane_flux( h,e );
        jn(:) = charge_neutral_flux(j, z );
       
        dcdt_v = intra_reaction_rate(c, s, k, z, v);
        
        for i = 1:n_comp
            dcdt_trans(:,i) = sum(jn(:,:,i), 2);
        end
        
        dcdt_v = dcdt_v - reshape(dcdt_trans, size(dcdt_v));
        
        dvdt_v  = dvdt(a, z, j, cap);
        
        state(1:(n_comp*n_species)) = reshape(dcdt_v, [1,(n_comp*n_species)]);
        state((n_comp*n_species)+1:end) = dvdt_v;
    end
end

