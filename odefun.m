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
%        size(n_comp, n_comp, n_species).
%    z: valence of species i, size(n_species)
%    o: volume of compartments, size(n_comp)
%    s: stochiometry matrix: (n_species, n_comp, n_reactions)
%    mu_ref: reference chemical potentail for species i in compartment a, 
%        size(n_comp, n_species)
%    k: forward and backward rates, n_reaction x 2.


    f = @fun;
    
    n_species = size(s, 1);
    n_comp = size(s, 2);
    n_reactions = size(s, 3);
    
    s = reshape(s, [n_species*n_comp, n_reactions]);
    
    % total membrane flus of species i from a to a'
    j = zeros(n_comp, n_comp, n_species);
    
    % charge-neutral flux of species from a' to a, 
    jt = zeros(n_comp, n_comp, n_species);
    

    for i=1:n_comp
        j(i,i,:) = 0;
    end
    
    % the stochiometric weighted sum of the reference chemical potentials. 
    delta_g = reference_chem_potential(k);
    
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
        
        %e = equilibrium_factor(c,z,v,l);
    
        %disp('z: ')
        %disp(z)
        %disp('e: ')
        %disp(e)

        %j(:,:,:) = membrane_flux(h,e);

        %disp('j: ')
        %disp(j)

        %jt(:,:,:) = charge_neutral_flux(j, z );

        %disp(jt)

        mu = electrochem_potential(c, z, v);

        % affinity, size(n_comp, n_reactions)
        y = affinity(s, delta_g, mu);

        % rate of reaction k in compartment a, size(n_comp, n_reactions)
        w = reaction_rate(r, y);

        dcdt_v = dcdt(o, a, jt, s, w);

        dvdt_v  = dvdt(a, z, j, cap);
        
        s(1:(n_comp*n_species)) = reshape(dcdt_v, [1,(n_comp*n_species)]);
        s((n_comp*n_species)+1:end) = dvdt_v;
    end
end

