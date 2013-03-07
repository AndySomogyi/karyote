function [ fun,  init ] = test(n_comp, n_species, n_reactions)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


%    cap: compartment capacitance, size(n_comp,n_comp), 
%        must have zero diagonal. 
cap = zero_diagonal(rand(n_comp, n_comp));

%    a: area between compartments (membrane), size(n_comp, n_comp), 
%        must have zero diagonal. 
a = zero_diagonal(rand(n_comp, n_comp));

%    l: thickness of membrane separating a and a' compartments, 
%        size(n_comp, n_comp), must have zero diagonal. 
l = zero_diagonal(rand(n_comp, n_comp));

%    h: membrane permeability for species i between compartments, 
%        size(n_comp, n_comp, n_species).
h = rand(n_comp, n_comp, n_species);

%    z: valence of species i, size(n_species)
z = randi([-3,3],[1,n_species]);

%    o: volume of compartments, size(n_comp)
o = rand(1,n_comp);

%    nu: = stoichiometric coefficient of species i for reaction k in
%        compartment a. dims: nreactions x ncomp x nspecies
nu = rand(n_reactions, n_comp, n_species);

%    mu_ref: reference chemical potentail for species i in compartment a, 
%        size(n_comp, n_species)
mu_ref = rand(n_comp, n_species);

%    r: reverse rate, size(n_comp x n_reactions)
r = rand(n_comp, n_reactions);

fun = odefun(cap, a, l, h, z, o, nu, mu_ref, r);

init = rand((n_comp * n_species) + n_comp, 1);


end

