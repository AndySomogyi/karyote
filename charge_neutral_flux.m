function [ jt ] = charge_neutral_flux(j, z )
% j: total flux, n_comp x n_comp x n_species
% z: is the vector of valences of the species that can transfer between 
%    compartment a1 and a2, dimensions: n_comp x n_comp x n_species

    sz = size(j);
    disp(sz)

    jt = zeros(sz);

    for i1=1:sz(1)
       for i2=1:sz(2)
           jt(i1,i2,:) =  squeeze(j(i1,i2,:))' - squeeze(z(i1,i2,:))' .* ...
               (...
               (squeeze(z(i1,i2,:))' * squeeze(j(i1,i2,:))) / ...
               (squeeze(z(i1,i2,:))' * squeeze(z(i1,i2,:))) ...
               );
       end
    end
    
    % verify that z . j == 0
    
    for i1=1:sz(1)
       for i2=1:sz(2)
           disp(squeeze(z(i1,i2,:))' * squeeze(jt(i1,i2,:)))
       end
    end
end

