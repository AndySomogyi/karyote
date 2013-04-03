function [mu] = electrochem_potential(c, z, v)
% Calculates the 'active' part of the chemical potential, 
% i.e. the chemical potential WITHOUT the reference chem potential. 
% The ref chem potential is added in later. 
%
% The membrane flux is divided into passive and deviatiric contributions
% mu_ref: reference chemical potential of species i in compartments, n compartments x m species
% c: concentrations, n comp x m species
% z: valences, m speies
% v: compartment voltages, n comp
% returns: mu, n_comp x n_species

    global R F T;
    mu = R * T *  log(c); % + F * v' * z;

end


