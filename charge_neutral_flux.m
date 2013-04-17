function [ jt ] = charge_neutral_flux(j, z_mobi)
% j: total flux, n_species x n_comp x n_comp
% z: is the vector of valences, dims: 1 x n_species

    sz = size(j);
    jt = zeros(sz);

    for i1=1:sz(2)
        for i2=1:sz(3)
            z = z_mobi(:,i1,i2)';
            zz = z * z';
            jj = j(:,i1,i2)';
            if zz ~= 0
                jt(:,i1,i2) = jj - z .* (z * jj')/(zz); 
            else
                jt(:,i1,i2) = jj;
            end           
        end
    end
end

