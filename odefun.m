function [ f ] = odefun(cap, a, l, h, z, o, s, k, r)
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
%    r: resitivity (for ionic current) between compartments. size(n_comp,n_comp).
%       the default value should be inf, and this needs to be a symmetric matrix, 
%       with inf diagonal. 


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
    assert(ismatrix(o) && length(o) == n_comp, 'length of compartment volume o must be n_comp');
    
    assert(length(o) ==  n_comp && ismatrix(o), 'volume vector must have the same length as n_comp');
    assert(ismatrix(cap) && size(cap,1) == n_comp && size(cap,1) == size(cap,2), ...
        'capacitance must be a square matrix with sides = n_comp');
    
    
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
    
    % calculate the inverse of the capacitance matrix. 
    % the given cap matrix is actually capacitance per unit area, so
    % this needs to be multiplied by membrane contact area to get 
    % the total membrane <-> membrane capacitance.
    % capacitance current convention is all currents leaving

    C = a .* cap;
    dvdt_inv = C(2:end,2:end);
    for i = 1:size(dvdt_inv, 1)
        dvdt_inv(i,i) = -sum(C(i+1,:));
    end
    clear C;
    dvdt_inv = inv(dvdt_inv);
    assert(all(all(isfinite(dvdt_inv))), 'capacitance matrix is not invertable');
    
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
        
        %ci = find(c <= 0);
        %if ~isempty(ci)   
        %    comp = floor((ci-1)/n_species)+1;
        %    species = mod(ci-1,n_species)+1;
        %    items = [comp',species'];
        %    warning('zero or negative concentration of [comp,species], [%s]', num2str(items));
        %end
        
        e = equilibrium_factor(c,z,v,l);
        j(:) = membrane_flux( h,e );
        jn(:) = charge_neutral_flux(j, z_mobi);
       
        dcdt_v = intra_reaction_rate(c, s, k, z, v);
        
        %fprintf('dcdt intra: %d \n', dcdt_v(273));
        
        for kk = 1:n_comp
            dcdt_trans(:,kk) = sum(jn(:,:,kk), 2) / o(kk);
        end
                
        dcdt_v = dcdt_v + reshape(dcdt_trans, size(dcdt_v));
        
        %fprintf('dcdt both: %d \n', dcdt_v(273));
        
        dvdt_v  = dvdt(a, z, j, dvdt_inv, v, r);
        %dvdt_v = zeros(size(v));
        %disp(dvdt_v);
        
        state(1:(n_comp*n_species)) = reshape(dcdt_v, [1,(n_comp*n_species)]);
        state((n_comp*n_species)+1:end) = dvdt_v;
    end
end

