function [e] = equilibrium_factor(c,z,v,l)
% c: concentration, n compartment x m species
% z: valences, m species
% v: voltages for compartments, n comp
% l: membrane thickness separating compartments, n comp x n comp
% = Faraday's constant Valence of i
% Electrical potential in Thickness of the membrane

    n_comp = length(v);
    e = zeros(n_comp, n_comp, length(z));
    nzi = logical(z ~= 0);
    zi = logical(z == 0);

    for a1 = 1:n_comp
        for a2 = 1:n_comp
            if a1 ~= a2
                % deal with zero and non zero valences
                e(a1,a2,nzi) = non_zero_z(...
                    c(a1,nzi),c(a2,nzi),z(nzi),v(a1)-v(a2),l(a1,a2));
                e(a1,a2,zi) = zero_z(c(a1,zi), c(a2,zi),l(a1,a2));
            end
        end
    end
end



function [e] = non_zero_z(c1, c2, z, dv, l)
    global F R T
    q = F * z * (dv) / (l * R * T);
    eql = exp(q * l);            
    e = q .* (                                ...
                  (c1 - c2 .* eql ) ...
                  ./                          ...
                  (1-eql)                     ...
             );
end

function [e] = zero_z(c1, c2, l)
    e = (c2 - c1)/l;
end

