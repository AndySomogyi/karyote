

function [mu] = electrochem_potential(mu_ref, c, z, v)
% The membrane flux is divided into passive and deviatiric contributions
% mu_ref: reference chemical potential of species i in compartments, n compartments x m species
% c: concentrations, n comp x m species
% z: valences, m speies
% v: compartment voltages, n comp

    mu = mu_ref + R * R * log(c) + F * z * v

end


