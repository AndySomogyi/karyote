function [e] = equilibrium_factor(c,z,v,l)
% c: concentration, size: (n compartment x m species), 1
% z: valences, m species
% v: voltages for compartments, n comp
% l: membrane thickness separating compartments, n comp x n comp
% = Faraday's constant Valence of i
% Electrical potential in Thickness of the membrane

    global F R T;

    n_comp = length(v);
    n_species = length(z);
    c = reshape(c, [n_species, n_comp]); % matlab uses colum major order
    e = zeros(length(z), n_comp, n_comp);
    

    for a1 = 1:n_comp
        for a2 = 1:n_comp
            if a1 ~= a2
                % deal with zero and non zero valences
                ll = l(a1,a2);
                dv = v(a1)-v(a2);
                q = F * z * (dv) / (ll * R * T);
                nzi = logical(q ~= 0);
                zi = logical(q == 0);
                e(nzi, a1,a2) = non_zero_q(q(nzi), c(nzi, a1),c(nzi, a2),ll);
                e(zi, a1,a2) = zero_q(c(zi, a1), c(zi, a2),ll);
            end
        end
    end
end



function [e] = non_zero_q(q, c1, c2, l)
    eql = exp(q * l);            
    e = q .* (                                ...
                  ((c1 - c2)' .* eql ) ...
                  ./                          ...
                  (1-eql)                     ...
             );
end

function [e] = zero_q(c1, c2, l)
    e = (c2 - c1)/l;
end

