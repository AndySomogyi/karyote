function [ c ] = neutralize_charge( c, z )
% concentration, n_species, n_comp
% z: 1,n_species;
    zp = zeros(size(z));
    lp = logical(z > 0);
    zp(lp) = z(lp);
    zm = zeros(size(z));
    lm = logical(z < 0);
    zm(lm) = z(lm);

    for i = 1:size(c,2)
        q = z * c(:,i);
        if q < 0
            scale = -q/(zp * c(:,i));
            c(lp,i) = (1+scale) * c(lp,i);
            warning('charge in compartment %i is %f, scaling positive charge by %f%', i, q, scale);
        elseif q > 0
            scale = -q/(zm * c(:,i));
            c(lm,i) = (1+scale) * c(lm,i);
            warning('charge in compartment %i is %f, scaling negative charge by %f%', i, q, scale);
        else
            fprintf('charge in compartment %i is 0\n', i);
        end 
    end
end

