function [ o ] = test()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    n_comp = 5;
    n_species = 10;
    
    c = zero_diagonal(rand(n_comp, n_species));
    v = rand(n_comp);
    a = zero_diagonal(rand(n_comp, n_comp));
    l = zero_diagonal(rand(n_comp, n_comp));
    nu = rand(n_comp, n_comp, n_species);
    h = rand(n_comp, n_comp, n_species);
    z = randi([-3,3],[1,n_species]);
    z_mob = randi([-3,3],[n_comp, n_comp ,n_species]);
    for i=1:n_comp
        z_mob(i,i,:) = 0;
    end
    
    j = rand(n_comp, n_comp ,n_species);
    for i=1:n_comp
        j(i,i,:) = 0;
    end
    




    e = equilibrium_factor(c,z,v,l);
    
    disp('z: ')
    disp(z)
    disp('e: ')
    disp(e)

    j = membrane_flux( h,e );
    
    
    
    disp('j: ')
    disp(j)

    jt = charge_neutral_flux(j, z_mob );
    
    disp(jt)




end

