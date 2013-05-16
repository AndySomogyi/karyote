function [ j ] = charge_neutral_flux(j, z)
% j: total flux, n_species x n_comp x n_comp
% z: is the vector of valences, dims: 1 x n_species

% concentration, n_species, n_comp
% z: 1,n_species;
    n_comp = length(z);
    zp = zeros(size(z));
    lp = logical(z > 0);
    zp(lp) = z(lp);
    zn = zeros(size(z));
    ln = logical(z < 0);
    zn(ln) = z(ln);
    
    for i = 1:n_comp
        for k=1:n_comp
            qp =  zp * abs(j(:,i,k));
            qn = -zn * abs(j(:,i,k));
                        
            if qp > qn
                s = qn / qp;
                j(lp,i,k) = s * j(lp,i,k);
            elseif qn > qp
                s = qp / qn;
                j(ln,i,k) = s * j(ln,i,k);
            end
            
            if z * abs(j(:,i,k)) > 1e-25
                qp =  zp * abs(j(:,i,k));
                qn = -zn * abs(j(:,i,k));

                if qp > qn
                    s = qn / qp;
                    j(lp,i,k) = s * j(lp,i,k);
                elseif qn > qp
                    s = qp / qn;
                    j(ln,i,k) = s * j(ln,i,k);
                end
                if z * abs(j(:,i,k)) > 1e-25
                    error('oh the horror, charge neutrality not reached, q: %d', z * abs(j(:,i,k)));
                end
            end
        end
    end
end

