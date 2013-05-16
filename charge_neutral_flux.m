function [ j ] = charge_neutral_flux(j, z)
% j: total flux, n_species x n_comp x n_comp
% z: is the vector of valences, dims: 1 x n_species

% concentration, n_species, n_comp
% z: 1,n_species;
    n_comp = size(j,2);
    zp = zeros(size(z));
    zn = zeros(size(z));
        
    for i = 1:n_comp
        for k=1:n_comp
            jj=j(:,i,k);
            ip1 = logical(jj > 0 & z' > 0); 
            ip2 = logical(jj < 0 & z' < 0);
            ip = ip1 | ip2;
            in1 = logical(jj > 0 & z' < 0);
            in2 = logical(jj < 0 & z' > 0);
            in = in1 | in2;
            
            zp(:) = 0;
            zp(ip1) =  z(ip1);
            zp(ip2) =  z(ip2);
            
            zn(:) = 0;
            zn(in1) = -z(in1);
            zn(in2) = -z(in2);
            
            qp =  zp * jj;
            qn =  zn * jj;
                        
            if qp > qn
                s = qn / qp;
                j(ip,i,k) = s * jj(ip);
            elseif qn > qp
                s = qp / qn;
                j(in,i,k) = s * jj(in);
            end
            
            q = z * j(:,i,k);
            
            if abs(q) > 1e-15
                error('death and horror!, charge neutral flux not reached, q:%d\n', q);
            end
          
        end
    end
end

