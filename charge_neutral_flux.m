function [ jt ] = charge_neutral_flux(j, z_mob )
% j: total flux, n_comp x n_comp x n_species
% z: is the vector of valences of the species that can transfer between 
%    compartment a1 and a2, dimensions: n_comp x n_comp x n_species

    sz = size(j);
    disp(sz)

    jt = zeros(sz);

    for i1=1:sz(1)
       for i2=1:sz(2)
           zz = squeeze(z_mob(i1,i2,:))' * squeeze(z_mob(i1,i2,:));
           if zz ~= 0
               jt(i1,i2,:) =  squeeze(j(i1,i2,:))' - squeeze(z_mob(i1,i2,:))' .* ...
                   (...
                   (squeeze(z_mob(i1,i2,:))' * squeeze(j(i1,i2,:))) / ...
                   (zz) ...
                   );
           else
               jt(i1,i2,:) = 0;
           end
       end
    end
    
    disp('jt: ')
    disp(jt)
    
    % verify that z . j == 0
    
    for i1=1:sz(1)
       for i2=1:sz(2)
           disp(squeeze(z_mob(i1,i2,:))' * squeeze(jt(i1,i2,:)))
       end
    end
end

