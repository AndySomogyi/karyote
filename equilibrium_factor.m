function [e] = equilibrium_factor(c,z,v,l)
% c: concentration, n compartment x m species
% z: valences, m species
% v: voltages for compartments, n comp
% l: membrane thickness separating compartments, n comp x n comp
% = Faraday's constant Valence of i
% Electrical potential in Thickness of the membrane

global F R T

    n_comp = length(v);
    e = zeros(n_comp, n_comp, length(z));

    for a1 = 1:n_comp
        for a2 = 1:n_comp
            if a1 ~= a2
                q = F * R * T * z * (v(a1) - v(a2)) / l(a1,a2);
                eql = exp(q * l(a1,a2));            
                e(a1,a2,:) = q .* (             ...
                    (c(a1,:) - c(a2,:) .* eql ) ...
                    ./                          ...
                    (1-eql)                     ...
                );
            end
        end
    end
end

