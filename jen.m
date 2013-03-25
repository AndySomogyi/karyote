% set the global constants
global F R T
F= 96485.3365; R=8.3144621; T=300;

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

%    initial values:
initc = zeros(n_comp, n_species);

%    initial voltage values:
initv = zeros(1,n_comp);


% now we can start assigning actual values


% capacitance is easy, they are all the same values, 
cap(:) = 2e-04;

% area between compartments, only compartments 1,2 are in contact, and
% comps 2 and 3 are in contact. 
a(1, 2) = 1.866e-4;
a(2, 3) = 1.866e-4;

% make symmetric
a = a' + triu(a,1);

% thickness of membrane in nm
%l(0,1) = 4;
l(1,2) = 2;
l(2,3) = 2;
% make symmetric
l = l' + triu(l,1);

% permeability, (n_comp, n_comp, n_species)
h(1,2,4) = 1;
h(2,3,4) = 1;

% valence of species
z(:) = [0 2 4 2 0 2 -1 -1 0];

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

nu(3,2,3) = -1;
nu(3,2,2) =  1;
nu(3,3,4) =  1;

nu(4,2,2) = -1;
nu(4,2,1) =  1;
nu(4,3,4) =  1;

nu(5,3,5) = -1;
nu(5,3,4) = -1;
nu(5,3,6) =  1;

nu(6,3,6) = -1;
nu(6,3,7) = -1;
nu(6,3,5) =  1;
nu(6,3,8) =  1;
nu(6,3,9) =  1;
nu(6,3,4) =  1;

%    mu_ref: reference chemical potentail, size(n_comp, n_species)
% how do I make all of the chemical potentials zero?

% r: reverse rate, size(n_comp x n_reactions)
% we need to do (comp,comp,reaction) I think for transmembrance reactions
r(2,1) = 1.0e11;
r(2,2) = 1.0e11;
r(3,3) = 1.0e11;
r(3,4) = 1.0e11;
r(3,5) = 1.0e10;
r(3,6) = 4.4693e2;

initc(:,1) = 1.0e-4;
initc(:,2) = 6.25e-4;
initc(:,3) = 6.0e-4;
initc(:,4) = 1.15e-4;
initc(:,5) = 2.7e-5;
initc(:,6) = 1.161e-3;
initc(:,7) = 4.89e-3;
initc(:,8) = 4.863e-3;
initc(:,9) = 2.7e-5;

% set initial voltage
initv(:) = [0.0 0.07 0.07];

% packed initial values
init = karyote_pack(initc, initv);

fun = odefun(cap, a, l, h, z, o, nu, mu_ref, r);

[t,y] = ode45(fun,[0 0.1], init);



