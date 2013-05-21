function [ c ] = neutralize_charge( c, z )
% concentration, n_species, n_comp
% z: 1,n_species;


    n_comp = size(c,2);
    zp = zeros(size(z));
    zn = zeros(size(z));
        
    for i = 1:n_comp

        cc=c(:,i);
        ip1 = logical(cc > 0 & z' > 0);
        ip2 = logical(cc < 0 & z' < 0);
        ip = ip1 | ip2;
        in1 = logical(cc > 0 & z' < 0);
        in2 = logical(cc < 0 & z' > 0);
        in = in1 | in2;
        
        zp(:) = 0;
        zp(ip1) =  z(ip1);
        zp(ip2) =  z(ip2);
        
        zn(:) = 0;
        zn(in1) = -z(in1);
        zn(in2) = -z(in2);
        
        qp =  zp * cc;
        qn =  zn * cc;
        
        if qp > qn
            s = qp / qn;
            c(in,i) = s * cc(in);
            fprintf('excess positive chage in comp %i, increasing negative charge by %3.2f%%\n', ...
                i, (s-1)*100);
        elseif qn > qp
            s = qn / qp;
            c(ip,i) = s * cc(ip);
            fprintf('excess negative chage in comp %i, increasing positive charge by %3.2f%%\n', ...
                i, (s-1)*100);
        else
            fprintf('charge is neutral in comp %i\n', i);
        end
        
        % the total charge caried by this flux, calculated with the
        % real valence vector, this had better be zero, or very near.
        q = z * c(:,i);
        
        fprintf('total charge in comp %i after neutralization: %d\n', i, q);
        
        if abs(q) > 1e-15
            error('death and horror!, unable to neutralize change %d in comp %i\n', ...
                q, i);
        end
    end
end

