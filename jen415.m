% set the global constants
global F R T
F= 96485.3365; R=8.3144621; T=300;

% define the size of the system
% the semicolon at the end is not required, all it does is suppress 
% diplaying of the value..
n_comp = 4;
n_species = 115;
n_reactions = 74;

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
h = zeros(n_species, n_comp, n_comp);

%    z: valence of species i, size(n_species)
z = zeros(1,n_species);

%    o: volume of compartments, size(n_comp)
o = zeros(1,n_comp);

%    nu: = stoichiometric coefficient of species i for reaction k in
%        compartment a. dims: nspecies, n_species x ncomp x n_reactions
nu = zeros(n_species, n_comp, n_reactions);

%    mu_ref: reference chemical potentail for species i in compartment a, 
%        size(n_comp, n_species)
mu_ref = zeros(n_comp, n_species);

%    k: rate (n_reaction, rate direction), 1 = forward, 2 = reverse
k = zeros(n_reactions, 2);

%    concentration inside c(n_comp x n_specie)
c = zeros(n_comp, n_species);

%    voltage
v = zeros(1,n_comp);

% "atom" count per species, not actually the real atom count, just
% how many things have to be conserved, for example, Adenesine is a 1
% because is always remains intact, but combines with other things.
ac = zeros(n_species);


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
% default is infinitly thick membrane
l(:) = Inf;
l(4,1) = 4;
l(1,2) = 2;
l(2,3) = 2;
% make symmetric
l = l' + triu(l,1);

% permeability, (n_comp, n_comp, n_species)
h(18,  4,1) = 1;
h(104, 4,1) = 1;
h(17,  1,2) = 1;
h(18,  1,2) = 1;
h(30,  1,2) = 1;
h(24,  1,2) = 1;
h(23,  1,2) = 1;
h(27,  1,2) = 1;
h(104, 1,2) = 1;
h(30,  2,1) = 1;
h(24,  2,1) = 1;
h(23,  2,1) = 1;
h(27,  2,1) = 1;
h(104, 2,1) = 1;
h(7,   2,3) = 1;
h(8,   2,3) = 1;
h(18,  2,3) = 1;
h(14,  2,3) = 1;
h(15,  2,3) = 1;
h(17,  2,3) = 1;
h(23,  2,3) = 1;
h(24,  2,3) = 1;
h(26,  2,3) = 1;
h(27,  2,3) = 1;
h(31,  2,3) = 1;
h(104, 2,3) = 1;
h(7,   3,2) = 1;
h(8,   3,2) = 1;
h(14,  3,2) = 1;
h(15,  3,2) = 1;
h(17,  3,2) = 1;
h(23,  3,2) = 1;
h(24,  3,2) = 1;
h(26,  3,2) = 1;
h(27,  3,2) = 1;
h(28,  3,2) = 1;
h(31,  3,2) = 1;
h(104, 3,2) = 1;

%h(:) = 0;

% valence of species
z(2) = -3;
z(3) = -3;
z(4) = -3;
z(5) = -2;
z(6) = -2;
z(7) = -2;
z(8) = -2;
z(9) = -2;
z(10) = -2;
z(14) =  1;
z(16) = -1;
z(17) =  1;
z(18) = -1;
z(21) =  2;
z(22) =  3;
z(23) = -4;
z(24) = -3;
z(25) = -2;
z(26) = -1;
z(27) =  1;
z(30) = -2;
z(31) = -1;
z(32) = -2;
z(41) = -1;
z(42) = -2;
z(43) = -4;
z(44) = -4;
z(46) = -2;
z(47) = -2;
z(48) = -3;
z(50) = -3;
z(51) = -3;
z(53) =  1;
z(55) =  4;
z(59) = -2;
z(64) = -2;
z(65) = -3;
z(66) = -1;
z(68) = -2;
z(69) = -2;
z(71) = -2;
z(72) = -2;
z(74) =  1;
z(75) = -1;
z(76) = -1;
z(77) =  2;
z(78) =  1;
z(80) =  1;
z(82) = -1;
z(84) = -2;
z(85) = -4;
z(86) = -4;
z(88) = -7;
z(92) =  1;
z(96) =  4;
z(98) = -1;
z(100) = -1;
z(101) = -1;
z(102) = -1;
z(103) = -1;
z(104) =  2;
z(105) =  2;
z(107) =  2;
z(108) = -1;
z(110) =  2;
z(111) =  4;
z(113) =  5;
z(114) =  1;
z(115) = -1;

% volume of each compartment
o(:) = [2.0e-7 7.0e07 1.6e-04, Inf];

% stoichometry, indexing: nspecies x ncomp x nreactions
% 1-34 are fast rxns, 35-44 are slow rxns (fast and slow are intracomparmental) and 45-74 are transmembrane rxns
nu(106,3,1) = -1;
nu(104,3,1) = -1;
nu(107,3,1) =  1;
nu(34,3,2) = -1;
nu(35,3,2) = -1;
nu(33,3,2) =  1;
nu(36,3,2) =  1;
nu(36,3,3) = -1;
nu(11,3,3) = -1;
nu(37,3,3) =  1;
nu(1,3,3) =  1;
nu(37,3,4) = -1;
nu(39,3,4) = -1;
nu(35,3,4) =  1;
nu(38,3,4) =  1;
nu(38,3,5) = -1;
nu(14,3,5) = -1;
nu(39,3,5) =  1;
nu(15,3,5) =  1;
nu(27,3,5) =  1;
nu(40,3,6) = -1;
nu(31,3,6) = -1;
nu(41,3,6) =  1;
%nu(41,3,7) = -1;
%nu(24,3,7) = -1;
%nu(43,3,7) =  1;
nu(44,3,8) = -1;
nu(42,3,8) =  1;
nu(25,3,8) =  1;
nu(42,3,9) = -1;
nu(40,3,9) =  1;
nu(30,3,9) =  1;
nu(45,3,10) = -1;
nu(10,3,10) = -1;
nu(46,3,10) =  1;
nu(46,3,11) = -1;
nu(1,3,11) = -1;
nu(47,3,11) = 1;
nu(48,3,12) = -1;
nu(45,3,12) = 1;
nu(2,3,12) = 1;
nu(49,3,13) = -1;
nu(2,3,13) = -1;
nu(50,3,13) = 1;
nu(51,3,14) = -1;
nu(29,3,14) = -1;
nu(49,3,14) = 1;
nu(4,3,14) = 1;
nu(52,3,15) = -1;
nu(14,3,15) = -1;
nu(53,3,15) =  1;
nu(53,3,16) = -1;
nu(102,3,16) = -1;
nu(54,3,16) = 1;
nu(55,3,17) = -1;
nu(52,3,17) = 1;
nu(15,3,17) = 1;
nu(104,3,17) = 2;
nu(56,3,18) = -1;
nu(104,3,18) = -1;
nu(105,3,18) = 1;
nu(57,3,19) = -1;
nu(58,3,19) = -1;
nu(105,3,19) = 1;
nu(59,3,19) = 1;
nu(59,3,20) = -1;
nu(11,3,20) = -1;
nu(60,3,20) = 1;
nu(6,3,20) = 1;
nu(60,3,21) = -1;
nu(62,3,21) = -1;
nu(58,3,21) = 1;
nu(61,3,21) = 1;
nu(61,3,22) = -1;
nu(14,3,22) = -1;
nu(62,3,22) = 1;
nu(15,3,22) = 1;
nu(27,3,22) = 1;
nu(63,3,23) = -1;
nu(6,3,23) = -1;
nu(64,3,23) = 1;
nu(65,3,24) = -1;
nu(66,3,24) = 1;
nu(7,3,24) = 1;
nu(66,3,25) = -1;
nu(24,3,25) = -1;
nu(63,3,25) = 1;
nu(23,3,25) = 1;
nu(70,3,26) = -1;
nu(8,3,26) = -1;
nu(71,3,26) = 1;
nu(72,3,27) = -1;
nu(70,3,27) = 1;
nu(9,3,27) = 1;
nu(73,3,28) = -1;
nu(14,3,28) = -1;
nu(74,3,28) = 1;
nu(74,3,29) = -1;
nu(9,3,29) = -1;
nu(75,3,29) = 1;
nu(76,3,30) = -1;
nu(78,3,30) = 1;
nu(10,3,30) = 1;
nu(78,3,31) = -1;
nu(73,3,31) = 1;
nu(15,3,31) = 1;
nu(27,3,31) = 1;
nu(19,3,32) = -1;
nu(18,3,32) = 1;
nu(27,3,32) = 1;
nu(23,3,33) = -1;
nu(77,3,33) = -1;
nu(30,3,33) = 1;
nu(24,3,34) = -1;
nu(77,3,34) = -1;
nu(31,3,34) = 1;
nu(107,3,35) = -1;
nu(108,3,35) = -1;
nu(106,3,35) = 1;
nu(26,3,35) = 1;
nu(33,3,35) = 1;
nu(104,3,35) = 1;
nu(18,3,36) = -1;
nu(33,3,36) = -1;
nu(20,3,36) = 1;
nu(34,3,36) = 1;
nu(103,3,36) = 1;
nu(43,3,37) = -1;
nu(44,3,37) = 1;
nu(47,3,38) = -1;
nu(48,3,38) = 1;
nu(101,3,38) = 1;
nu(27,3,38) = 2;
nu(50,3,39) = -1;
nu(51,3,39) = 1;
nu(29,3,39) = 1;
nu(54,3,40) = -1;
nu(104,3,40) = -2;
nu(55,3,40) = 1;
nu(5,3,40) = 1;
nu(32,3,40) = 1;
nu(27,3,40) = 2;
nu(77,3,40) = 1;
nu(5,3,41) = -1;
nu(105,3,41) = -1;
nu(57,3,41) = 1;
nu(20,3,41) = 1;
nu(64,3,42) = -1;
nu(26,3,42) = -1;
nu(65,3,42) = 1;
nu(11,3,42) = 1;
nu(71,3,43) = -1;
nu(29,3,43) = -1;
nu(72,3,43) = 1;
nu(75,3,44) = -1;
nu(76,3,44) = 1;
nu(27,1,45) = -1;
nu(79,2,45) = -1;
nu(80,2,45) = 1;
nu(80,2,46) = -1;
nu(18,1,46) = -1;
nu(81,2,46) = 1;
nu(81,2,47) = -1;
nu(18,3,47) = -1;
nu(80,2,47) = -1;
nu(80,2,48) = -1;
nu(79,2,48) = 1;
nu(27,3,48) = 1;
nu(84,2,49) = -1;
nu(31,3,49) = -1;
nu(26,3,49) = -1;
nu(85,2,49) = 1;
nu(85,2,50) = -1;
nu(27,1,50) = -4;
nu(86,2,50) = 1;
nu(27,3,50) = 4;
nu(86,2,51) = -1;
nu(84,2,51) = 1;
nu(30,3,51) = 1;
nu(87,2,52) = -1;
nu(24,1,52) = -1;
nu(23,3,52) = -1;
nu(88,2,52) = 1;
nu(88,2,53) = -1;
nu(24,3,53) = 1;
nu(23,1,53) = 1;
nu(89,2,54) = -1;
nu(15,3,54) = -1;
nu(90,2,54) = 1;
nu(90,2,55) = -1;
nu(12,2,55) = -1;
nu(91,2,55) = 1;
nu(91,2,56) = -1;
nu(27,3,56) = -5;
nu(92,2,56) = 1;
nu(13,2,56) = 1;
nu(27,1,56) = 4;
nu(92,2,57) = -1;
nu(89,2,57) = 1;
nu(14,3,57) = 1;
nu(67,2,58) = -1;
nu(7,3,58) = -1;
nu(68,2,58) = 1;
nu(68,2,59) = -1;
nu(12,2,59) = -1;
nu(69,2,59) = 1;
nu(13,2,59) = 1;
nu(69,2,60) = -1;
nu(67,2,60) = 1;
nu(8,3,60) = 1;
nu(93,2,61) = -1;
nu(13,2,61) = -1;
nu(22,1,61) = -1;
nu(94,2,61) = 1;
nu(21,1,61) = 1;
nu(27,1,61) = 2;
nu(103,1,61) = 1;
nu(94,2,62) = -1;
nu(13,2,62) = -1;
nu(22,1,62) = -1;
nu(27,3,62) = -2;
nu(103,1,62) = -1;
nu(93,2,62) = 1;
nu(12,2,62) = 1;
nu(13,2,62) = 1;
nu(21,2,62) = 1;
nu(27,1,62) = 2;
nu(95,2,63) = -1;
nu(21,1,63) = -2;
nu(96,2,63) = 1;
nu(96,2,64) = -1;
nu(28,3,64) = -1;
nu(97,2,64) = 1;
nu(22,1,64) = 2;
nu(103,1,64) = 2;
nu(97,2,65) = -1;
nu(21,1,65) = -1;
nu(27,3,65) = -2;
nu(103,3,65) = -1;
nu(98,2,65) = 1;
nu(22,1,65) = 1;
nu(27,1,65) = 1;
nu(98,2,66) = -1;
nu(21,1,66) = -1;
nu(27,3,66) = -2;
nu(99,2,66) = 1;
nu(22,1,66) = 1;
nu(27,1,66) = 1;
nu(103,3,66) = 1;
nu(99,2,67) = -1;
nu(27,3,67) = -2;
nu(103,1,67) = -2;
nu(100,2,67) = 1;
nu(29,3,67) = 1;
nu(27,1,67) = 1;
nu(100,2,68) = -1;
nu(27,3,68) = -2;
nu(95,2,68) = 1;
nu(29,3,68) = 1;
nu(27,1,68) = 1;
nu(109,2,69) = -1;
nu(104,1,69) = -1;
nu(110,2,69) = 1;
nu(110,2,70) = -1;
nu(104,1,70) = -1;
nu(111,2,70) = 1;
nu(111,2,71) = -1;
nu(110,2,71) = 1;
nu(104,3,71) = 1;
nu(110,2,72) = -1;
nu(109,2,72) = 1;
nu(104,3,72) = 1;
nu(112,2,73) = -1;
nu(104,3,73) = -1;
nu(17,1,73) = -3;
nu(113,2,73) = 1;
nu(113,2,74) = -1;
nu(112,3,74) = 1;
nu(104,1,74) = 1;
nu(17,2,74) = 3;

% only do intra compartment reactions
nu = nu(:,:,1:44);

%     k = (n_reactions, 2);
k(1,1) = 1.0e1;
k(2,1) = 1.0e1;
k(3,1) = 1.0e1;
k(4,1) = 1.0e1;
k(5,1) = 1.0e1;
k(6,1) = 1.0e1;
k(7,1) = 1.0e1;
k(8,1) = 1.0e1;
k(9,1) = 1.0e1;
k(10,1) = 1.0e1;
k(11,1) = 1.0e1;
k(12,1) = 1.0e1;
k(13,1) = 1.0e1;
k(14,1) = 1.0e1;
k(15,1) = 1.0e1;
k(16,1) = 1.0e1;
k(17,1) = 1.0e1;
k(18,1) = 1.0e1;
k(19,1) = 1.0e1;
k(20,1) = 1.0e1;
k(21,1) = 1.0e1;
k(22,1) = 1.0e1;
k(23,1) = 1.0e1;
k(24,1) = 1.0e1;
k(25,1) = 1.0e1;
k(26,1) = 1.0e1;
k(27,1) = 1.0e1;
k(28,1) = 1.0e1;
k(29,1) = 1.0e1;
k(30,1) = 1.0e1;
k(31,1) = 1.0e1;
k(32,1) = 1.0e1;
k(33,1) = 1.0e1;
k(34,1) = 1.0e1;
k(35,1) = 4.4693e2;
k(36,1) = 4.4693e2;
k(37,1) = 6.7615e6;
k(38,1) = 1.15473e5;
k(39,1) = 3.474e3;
k(40,1) = 1.905286e6;
k(41,1) = 1.7314583e7;
k(42,1) = 3.266667e6;
k(43,1) = 9.7e5;
k(44,1) = 5.47857e5;
k(45,1) = 1.22e5;
k(46,1) = 1.0e11;
k(47,1) = 1.0e11;
k(48,1) = 1.22e5;
k(49,1) = 1.0e11;
k(50,1) = 6.889e5;
k(51,1) = 1.0e11;
k(52,1) = 1.0e11;
k(53,1) = 3.7453;
k(54,1) = 1.0e11;
k(55,1) = 1.0e11;
k(56,1) = 9.2083e6;
k(57,1) = 1.0e11;
k(58,1) = 5.567e5;
k(59,1) = 1.0e11;
k(60,1) = 1.0e11;
k(61,1) = 1.0e11;
k(62,1) = 1.03e6;
k(63,1) = 1.0e11;
k(64,1) = 6.0e9;
k(65,1) = 1.0e11;
k(66,1) = 1.0e11;
k(67,1) = 1.0e11;
k(68,1) = 1.0e11;
k(69,1) = 1.0e11;
k(70,1) = 1.0e11;
k(71,1) = 1.0e5;
k(72,1) = 1.0e11;
k(73,1) = 1.0e11;
k(74,1) = 1.0e5;
k(:,2) = 1.0e-2;

k = k(1:44,:);

% concentration inside c(n_comp x n_specie)
% this should actually be n_species, n_comp, but already writen
% as n_comp, n_species, so in order to save time, 
% we'll just take the transpose of it:)

% first set each concentration to a very small value, can not be zero.
c(:) = 1.0e-10;

c(1,16) = 9.6880398e-3;
c(1,17) = 3.40570052e-2;
c(1,18) = 1.87e-4;
c(1,21) = 2.05e-3;
c(1,22) = 2.73333333333333e-3;
c(1,23) = 3.535e-3;
c(1,24) = 1.7e-3;
c(1,25) = 1.51e-4;
c(1,27) = 3.98e-8;
c(1,30) = 7.07e-3;
c(1,31) = 5.1e-3;
c(1,103) = 5.2e-9;
c(1,104) = 1.15e-3;
c(1,114) = 5.20417042793042e-18;
c(2,12) = 1.9e-2;
c(2,13) = 1.0e-3;
c(2,67) = 8.9e-5;
c(2,68) = 1.11095e-3;
c(2,69) = 1.525e-4;
c(2,79) = 1.0e-4;
c(2,80) = 1.000052e-4;
c(2,81) = 2.870052e-4;
c(2,82) = 2.87e-4;
c(2,83) = 3.17e-4;
c(2,84) = 3.6935e-3;
c(2,85) = 3.4875e-3;
c(2,86) = 2.27175e-3;
c(2,87) = 1.0e-4;
c(2,88) = 7.96715788571429e-3;
c(2,89) = 5.29e-5;
c(2,90) = 1.0529e-3;
c(2,91) = 2.00529e-2;
c(2,92) = 2.0529e-3;
c(2,93) = 1.11e-4;
c(2,94) = 1.9111e-2;
c(2,95) = 3.17e-4;
c(2,96) = 1.62330237e-2;
c(2,97) = 8.517e-3;
c(2,98) = 8.517e-3;
c(2,99) = 8.517e-3;
c(2,100) = 8.517e-3;
c(2,109) = 1.0e-4;
c(2,110) = 6.25e-4;
c(2,111) = 6.0e-4;
c(2,112) = 1.0e-4;
c(2,113) = 7.06140104e-3;
c(3,1) = 3.9e-6;
c(3,2) = 3.75e-4;
c(3,3) = 3.75e-4;
c(3,4) = 2.9e-5;
c(3,5) = 2.02e-4;
c(3,6) = 6.3e-4;
c(3,7) = 1.065e-3;
c(3,8) = 1.08e-4;
c(3,9) = 4.91e-4;
c(3,10) = 1.0e-5;
c(3,11) = 1.8e-4;
c(3,14) = 2.0e-3;
c(3,15) = 1.0e-3;
c(3,16) = 1.014060052e-1;
c(3,17) = 1.860538e-1;
c(3,18) = 1.87e-4;
c(3,19) = 1.87e-4;
c(3,20) = 2.14e-4;
c(3,23) = 3.535e-3;
c(3,24) = 1.7e-3;
c(3,25) = 1.51e-4;
c(3,26) = 4.863e-3;
c(3,27) = 5.2e-9;
c(3,28) = 1.0e-6;
c(3,29) = 6.435e-1;
c(3,30) = 3.535e-3;
c(3,31) = 1.7e-3;
c(3,32) = 1.07e-4;
c(3,33) = 2.7e-5;
c(3,34) = 2.7e-5;
c(3,35) = 2.7e-5;
c(3,36) = 4.28e-4;
c(3,37) = 2.7e-5;
c(3,38) = 2.7e-5;
c(3,40) = 1.0e-4;
c(3,41) = 7.17e-3;
c(3,42) = 9e-3;
c(3,43) = 1.725e-3;
c(3,44) = 1.868e-3;
c(3,45) = 1.1e-5;
c(3,46) = 1.05e-5;
c(3,47) = 1.245e-5;
c(3,48) = 3.78666666666667e-4;
c(3,49) = 1e-4;
c(3,50) = 4.08333333333333e-4;
c(3,51) = 4.08333333333333e-4;
c(3,52) = 4.7e-5;
c(3,53) = 2.047e-3;
c(3,54) = 2.134e-3;
c(3,55) = 1.047e-3;
c(3,56) = 2.2e-5;
c(3,57) = 2.322e-3;
c(3,58) = 2.2e-5;
c(3,59) = 1.076e-3;
c(3,60) = 2.2e-5;
c(3,61) = 2.2e-5;
c(3,62) = 2.2e-5;
c(3,63) = 1.08e-4;
c(3,64) = 6.84e-4;
c(3,65) = 2.077e-3;
c(3,66) = 4.971e-3;
c(3,70) = 5e-4;
c(3,71) = 3.58e-4;
c(3,72) = 7.41e-4;
c(3,73) = 9.1e-5;
c(3,74) = 2.091e-3;
c(3,75) = 3.073e-3;
c(3,76) = 1.1110052e-3;
c(3,77) = 1.0e-3;
c(3,78) = 1.0910052e-3;
c(3,101) = 1.8e-4;
c(3,102) = 2.087e-3;
c(3,104) = 1.15e-3;
c(3,105) = 1.161e-3;
c(3,106) = 2.7e-5;
c(3,107) = 1.161e-3;
c(3,108) = 4.89e-3;
c(3,114) = 1.202005498765e-2;
c(3,115) = 4.668955e-3;
c(4,18) = 2.3000052e-3;
c(4,19) = 1.87e-4;
c(4,27) = 5.20e-9;
c(4,104) = 1.15e-3;

% just take the tanspose of it...
c=c';


% voltage v(n_comp)
% initial values for compartment potentials. 
v(:) = 0;

% atom count
ac(1) = 
ac(2)


verify_stochiometry(nu, z);

% make the function that the integrator calls. 
fun = odefun(cap, a, l, h, z, o, nu, k);

% pack the initial values into the state vector
state = karyote_pack(c,v);

yp0 = zeros(size(state));

% test the function, call it once with the starting state vector. 
%options = odeset('NonNegative', 1:length(state)-n_comp, ...
%        'MaxStep', 1e-5);
    
% test the function, call it once with the starting state vector. 
t0 = 0;
tf = 0.00001;
options = odeset('NonNegative', 1:(length(state)-n_comp), ...
                 'RelTol', 1e-15, ...
                 'AbsTol', 1e-30, ...
                 'InitialStep', 0.01*abs(t0-tf));
[t,y] = ode15s(fun, [t0 tf], state, options);

c = y(:,2*n_species+1:3*n_species);

disp(c(1,:));
disp(c(end,:));

v = y(:,end-3:end);

plot(t,c);


