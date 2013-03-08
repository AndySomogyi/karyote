% define the size of the system
% the semicolon at the end is not required, all it does is suppress 
% diplaying of the value..
n_comp = 3;
n_species = 9;
n_reactions = 6;

% first make empty (zero) matricies to store the parameters, 
% easier this way as most parameters are zero.

%    cap: compartment capacitance, size(n_comp,n_comp), 
%        must have zero diagonal. 
cap = zeros(n_comp, n_comp);

%    a: area between compartments (membrane), size(n_comp, n_comp), 
%        must have zero diagonal. 
a = zeros(n_comp, n_comp);

%    l: thickness of membrane separating a and a' compartments, 
%        size(n_comp, n_comp), must have zero diagonal. 
l = zeros(n_comp, n_comp);

%    h: membrane permeability for species i between compartments, 
%        size(n_comp, n_comp, n_species).
h = zeros(n_comp, n_comp, n_species);

%    z: valence of species i, size(n_species)
z = zeros(1,n_species);

%    o: volume of compartments, size(n_comp)
o = zeros(1,n_comp);

%    nu: = stoichiometric coefficient of species i for reaction k in
%        compartment a. dims: nreactions x ncomp x nspecies
nu = zeros(n_reactions, n_comp, n_species);

%    mu_ref: reference chemical potentail for species i in compartment a, 
%        size(n_comp, n_species)
mu_ref = zeros(n_comp, n_species);

%    r: reverse rate, size(n_comp x n_reactions)
r = zeros(n_comp, n_reactions);


% now we can start assigning actual values


% capacitance is easy, they are all the same values, 
cap(:) = 2e-04;

% area between compartments, only compartments 1,2 are in contact, and
% comps 2 and 3 are in contact. 
a(1, 2) = 1.866e-4;
a(2, 1) = 1.866e-4;
a(2, 3) = 1.866e-4;
a(3, 2) = 1.866e-4;

% volume of each compartment
o(:) = [2.0e-7 7.0e07 1.6e-04];

% stoichometry, indexing: nreactions x ncomp x nspecies
% MCU(2) + Ca2+(1) <-> MCU[Ca2+](2)
nu(1,2,1) = -1;
nu(1,1,4) = -1;
nu(1,2,2) =  1;

nu(2,2,2) = -1;
nu(2,1,4) = -1;
nu(2,2,3) =  1;








