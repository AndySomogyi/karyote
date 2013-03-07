function [ f ] = odefun(cap, a, l, h, z, o, nu, mu_ref, r)
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
%    nu: = stoichiometric coefficient of species i for reaction k in
%        compartment a. dims: nreactions x ncomp x nspecies
%    mu_ref: reference chemical potentail for species i in compartment a, 
%        size(n_comp, n_species)
%    r: reverse rate, size(n_comp x n_speces)


    f = @fun;

    n_comp = size(cap, 1);
    n_species = length(z);
    
    % total membrane flus of species i from a to a'
    j = zeros(n_comp, n_comp, n_species);
    
    % charge-neutral flux of species from a' to a, 
    jt = zeros(n_comp, n_comp, n_species);
    

    for i=1:n_comp
        j(i,i,:) = 0;
    end
    
    function [s] = fun(t,s)
        % the entire state, [concentration; potential] is wrapped
        % up in the state vector s. s should be 
        % (n_comp*n_species) + n_comp long 
        
        % concentration
        % c: size(n_comp, n_comp)
        %
        % voltage
        % v: size(1,n_comp)
        
        fprintf('time: %d\n', t);
        
        c = reshape(s(1:(n_comp*n_species)), [n_comp, n_species]);
        v = reshape(s((n_comp*n_species)+1:end), [1, n_comp]);
        
        e = equilibrium_factor(c,z,v,l);
    
        %disp('z: ')
        %disp(z)
        %disp('e: ')
        %disp(e)

        j(:,:,:) = membrane_flux(h,e);

        %disp('j: ')
        %disp(j)

        jt(:,:,:) = charge_neutral_flux(j, z );

        %disp(jt)

        mu = electrochem_potential(mu_ref, c, z, v);

        % affinity, size(n_comp, n_reactions)
        y = affinity(nu, mu);

        % rate of reaction k in compartment a, size(n_comp, n_reactions)
        w = reaction_rate(r, y);

        dcdt_v = dcdt(o, a, jt, nu, w);

        dvdt_v  = dvdt(a, z, j, cap);
        
        s(1:(n_comp*n_species)) = reshape(dcdt_v, [1,(n_comp*n_species)]);
        s((n_comp*n_species)+1:end) = dvdt_v;
    end
end

