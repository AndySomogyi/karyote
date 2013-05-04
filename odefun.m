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
    
    fprintf('creating system for %i species, %i compartments, %i reactions\n', ...
        n_species, n_comp, n_reactions);
    
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
    
    assert(length(o) ==  n_comp && ismatrix(o), 'volume vector must have the same length as n_comp');
    
    % the mobile valences, these are species
    % which can either diffuse or are actively pumped between compartments.
    % start with species that can diffuse, put in active pumping later.
    z_mobi = zeros(size(h));
    for ii=1:n_comp
        for jj=1:n_comp
            nz = logical(h(:,ii,jj) > 0);
            z_mobi(nz) = z(nz);
        end
    end
    
    function [state] = fun(t,state)
        % the entire state, [concentration; potential] is wrapped
        % up in the state vector s. s should be 
        % (n_comp*n_species) + n_comp long 
        
        % concentration
        % c: size(n_comp, n_comp)
        %
        % voltage
        % v: size(1,n_comp)
        
        %fprintf('time: %d\n', t);
        
        % concentration vector, organized by compartment, then species. 
        c = state(1:(n_comp*n_species));
        % voltage vector, n_comp x 1
        v = state((n_comp*n_species)+1:end);
        
        % check to see that concentration is positive, and NOT 0 or
        % negative.
        
        %fprintf('value: %d\n', c(273))
        
        %c = abs(c);
        
        ci = find(c <= 0);
        if ~isempty(ci)   
            comp = floor((ci-1)/n_species)+1;
            species = mod(ci-1,n_species)+1;
            items = [comp',species'];
            error('zero or negative concentration of [comp,species], [%s]', num2str(items));
        end
        
        e = equilibrium_factor(c,z,v,l);
        j(:) = membrane_flux( h,e );
        jn(:) = charge_neutral_flux(j, z_mobi);
       
        dcdt_v = intra_reaction_rate(c, s, k, z, v);
        
        %fprintf('dcdt intra: %d \n', dcdt_v(273));
        
        for i = 1:n_comp
            dcdt_trans(:,i) = (1/o(i)) * sum(jn(:,:,i), 2);
        end
                
        dcdt_v = dcdt_v - reshape(dcdt_trans, size(dcdt_v));
        
        %fprintf('dcdt both: %d \n', dcdt_v(273));
        
        
        
        dvdt_v  = dvdt(a, z, j, cap);
        
        state(1:(n_comp*n_species)) = reshape(dcdt_v, [1,(n_comp*n_species)]);
        state((n_comp*n_species)+1:end) = dvdt_v;
    end
end

