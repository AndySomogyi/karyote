% set the global constants
global F R T
F= 96485.3365; R=8.3144621; T=300;

% define the size of the system
% the semicolon at the end is not required, all it does is suppress 
% diplaying of the value..
n_comp = 4;
n_species = 115;
n_intra_reactions = 44;
n_trans_reactions = 30;

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

%    intra-compartment stoichiometry matrix, 
%    stoichiometric coefficient of species i for reaction k in
%    compartment a. size: nspecies, n_species x ncomp x n_reactions
si = zeros(n_species, n_comp, n_intra_reactions);

%    trans-compartment stoichiometry matrix, 
%    stoichiometric coefficient of species i for reaction k in
%    compartment a. size: nspecies, n_species x ncomp x n_reactions
st = zeros(n_species, n_comp, n_trans_reactions);

%    mu_ref: reference chemical potentail for species i in compartment a, 
%        size(n_comp, n_species)
mu_ref = zeros(n_comp, n_species);

%    ki: intra-compartment reaction rates,
%    size: (n_reaction, rate direction), 1 = forward, 2 = reverse
ki = zeros(n_intra_reactions, 2);

%    kt: trans-compartment reaction rates,
%    size: (n_reaction, rate direction), 1 = forward, 2 = reverse
kt = zeros(n_trans_reactions, 2);

%    initial concentration, size(n_species, n_comp)
c0 = zeros(n_species, n_comp);

%    initial voltage values
v0 = zeros(1,n_comp);

%    resitivity between compartments, in most cases, should be inf
r = inf(n_comp, n_comp); 

% "atom" count per species, not actually the real atom count, just
% how many things have to be conserved, for example, Adenesine is a 1
% because is always remains intact, but combines with other things.
ac = zeros(n_species);


% now we can start assigning actual values


% capacitance is easy, they are all the same values, 
cap(:) = 2e-04;

% area between compartments, only compartments 1,2 are in contact, and
% comps 2 and 3 are in contact. 
a(1, 4) = 1.866e-4;
a(2, 3) = 1.866e-4;
a(3, 4) = 1.866e-4;

% make symmetric
a = a' + triu(a,1);

% thickness of membrane in nm
% default is infinitly thick membrane
l(:) = Inf;

l(4,1) = 4;
l(2,3) = 2;
l(3,4) = 2;

l(1,4) = 4;
l(3,2) = 2;
l(4,3) = 2;


%% permeability, (n_species, n_comp, n_comp) %%
% diffuse from universe to into outer comp (1)
h(18,  1,4) = 1;
h(104, 1,4) = 1;
h(17,  2,3) = 1;
h(18,  2,3) = 1;
h(30,  2,3) = 1;
h(24,  2,3) = 1;
h(23,  2,3) = 1;
h(27,  2,3) = 1;
h(104, 2,3) = 1;
h(30,  3,1) = 1;
h(24,  3,1) = 1;
h(23,  3,1) = 1;
h(27,  3,1) = 1;
h(104, 3,1) = 1;
h(7,   3,4) = 1;
h(8,   3,4) = 1;
h(18,  3,4) = 1;
h(14,  3,4) = 1;
h(15,  3,4) = 1;
h(17,  3,4) = 1;
h(23,  3,4) = 1;
h(24,  3,4) = 1;
h(26,  3,4) = 1;
h(27,  3,4) = 1;
h(31,  3,4) = 1;
h(104, 3,4) = 1;
h(7,   4,3) = 1;
h(8,   4,3) = 1;
h(14,  4,3) = 1;
h(15,  4,3) = 1;
h(17,  4,3) = 1;
h(23,  4,3) = 1;
h(24,  4,3) = 1;
h(26,  4,3) = 1;
h(27,  4,3) = 1;
h(28,  4,3) = 1;
h(31,  4,3) = 1;
h(104, 4,3) = 1;


%% valence of species %%
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

%% volume of each compartment
o(:) = [2.0e-7 7.0e07 1.6e-04, Inf];

%% stoichometry of intra-compartment reactions, indexing: nspecies x ncomp x nreactions
% 1-34 are fast rxns, 35-44 are slow rxns (fast and slow are intracomparmental) 
si(106,4,1) = -1;
si(104,4,1) = -1;
si(107,4,1) =  1;
si(34,4,2) = -1;
si(35,4,2) = -1;
si(33,4,2) =  1;
si(36,4,2) =  1;
si(36,4,3) = -1;
si(11,4,3) = -1;
si(37,4,3) =  1;
si(1,4,3) =  1;
si(37,4,4) = -1;
si(39,4,4) = -1;
si(35,4,4) =  1;
si(38,4,4) =  1;
si(38,4,5) = -1;
si(14,4,5) = -1;
si(39,4,5) =  1;
si(15,4,5) =  1;
si(27,4,5) =  1;
si(40,4,6) = -1;
si(31,4,6) = -1;
si(41,4,6) =  1;
si(41,4,7) = -1;
si(24,4,7) = -1;
si(43,4,7) =  1;
si(44,4,8) = -1;
si(42,4,8) =  1;
si(25,4,8) =  1;
si(42,4,9) = -1;
si(40,4,9) =  1;
si(30,4,9) =  1;
si(45,4,10) = -1;
si(10,4,10) = -1;
si(46,4,10) =  1;
si(46,4,11) = -1;
si(1,4,11) = -1;
si(47,4,11) = 1;
si(48,4,12) = -1;
si(45,4,12) = 1;
si(2,4,12) = 1;
si(49,4,13) = -1;
si(2,4,13) = -1;
si(50,4,13) = 1;
si(51,4,14) = -1;
si(29,4,14) = -1;
si(49,4,14) = 1;
si(4,4,14) = 1;
si(52,4,15) = -1;
si(14,4,15) = -1;
si(53,4,15) =  1;
si(53,4,16) = -1;
si(102,4,16) = -1;
si(54,4,16) = 1;
si(55,4,17) = -1;
si(52,4,17) = 1;
si(15,4,17) = 1;
si(104,4,17) = 2;
si(56,4,18) = -1;
si(104,4,18) = -1;
si(105,4,18) = 1;
si(57,4,19) = -1;
si(58,4,19) = -1;
si(105,4,19) = 1;
si(59,4,19) = 1;
si(59,4,20) = -1;
si(11,4,20) = -1;
si(60,4,20) = 1;
si(6,4,20) = 1;
si(60,4,21) = -1;
si(62,4,21) = -1;
si(58,4,21) = 1;
si(61,4,21) = 1;
si(61,4,22) = -1;
si(14,4,22) = -1;
si(62,4,22) = 1;
si(15,4,22) = 1;
si(27,4,22) = 1;
si(63,4,23) = -1;
si(6,4,23) = -1;
si(64,4,23) = 1;
si(65,4,24) = -1;
si(66,4,24) = 1;
si(7,4,24) = 1;
si(66,4,25) = -1;
si(24,4,25) = -1;
si(63,4,25) = 1;
si(23,4,25) = 1;
si(70,4,26) = -1;
si(8,4,26) = -1;
si(71,4,26) = 1;
si(72,4,27) = -1;
si(70,4,27) = 1;
si(9,4,27) = 1;
si(73,4,28) = -1;
si(14,4,28) = -1;
si(74,4,28) = 1;
si(74,4,29) = -1;
si(9,4,29) = -1;
si(75,4,29) = 1;
si(76,4,30) = -1;
si(78,4,30) = 1;
si(10,4,30) = 1;
si(78,4,31) = -1;
si(73,4,31) = 1;
si(15,4,31) = 1;
si(27,4,31) = 1;
si(19,4,32) = -1;
si(18,4,32) = 1;
si(27,4,32) = 1;
si(23,4,33) = -1;
si(77,4,33) = -1;
si(30,4,33) = 1;
si(24,4,34) = -1;
si(77,4,34) = -1;
si(31,4,34) = 1;
si(107,4,35) = -1;
si(108,4,35) = -1;
si(106,4,35) = 1;
si(26,4,35) = 1;
si(33,4,35) = 1;
si(104,4,35) = 1;
si(18,4,36) = -1;
si(33,4,36) = -1;
si(27,4,36) = -1;
si(20,4,36) = 1;
si(34,4,36) = 1;
si(43,4,37) = -1;
si(44,4,37) = 1;
si(47,4,38) = -1;
si(48,4,38) = 1;
si(101,4,38) = 1;
si(27,4,38) = 2;
si(50,4,39) = -1;
si(51,4,39) = 1;
si(29,4,39) = 1;
si(54,4,40) = -1;
si(104,4,40) = -2;
si(55,4,40) = 1;
si(5,4,40) = 1;
si(32,4,40) = 1;
si(27,4,40) = 2;
si(77,4,40) = 1;
si(5,4,41) = -1;
si(105,4,41) = -1;
si(57,4,41) = 1;
si(20,4,41) = 1;
si(64,4,42) = -1;
si(26,4,42) = -1;
si(65,4,42) = 1;
si(11,4,42) = 1;
si(71,4,43) = -1;
si(29,4,43) = -1;
si(72,4,43) = 1;
si(75,4,44) = -1;
si(76,4,44) = 1;


%% stoichiomery matrix for trans compartment reactions 
%% and 45-74 are transmembrane rxns
st(27,2,1) = -1;
st(79,3,1) = -1;
st(80,3,1) = 1;
st(80,3,2) = -1;
st(18,2,2) = -1;
st(81,3,2) = 1;
st(81,3,3) = -1;
st(18,4,3) = -1;
st(80,3,3) = -1;
st(80,3,4) = -1;
st(79,3,4) = 1;
st(27,4,4) = 1;
st(84,3,5) = -1;
st(31,4,5) = -1;
st(26,4,5) = -1;
st(85,3,5) = 1;
st(85,3,6) = -1;
st(27,2,6) = -4;
st(86,3,6) = 1;
st(27,4,6) = 4;
st(86,3,7) = -1;
st(84,3,7) = 1;
st(30,4,7) = 1;
st(87,3,8) = -1;
st(24,2,8) = -1;
st(23,4,8) = -1;
st(88,3,8) = 1;
st(88,3,9) = -1;
st(24,4,9) = 1;
st(23,2,9) = 1;
st(89,3,10) = -1;
st(15,4,10) = -1;
st(90,3,10) = 1;
st(90,3,11) = -1;
st(12,3,11) = -1;
st(91,3,11) = 1;
st(91,3,12) = -1;
st(27,4,12) = -5;
st(92,3,12) = 1;
st(13,3,12) = 1;
st(27,2,12) = 4;
st(92,3,13) = -1;
st(89,3,13) = 1;
st(14,4,13) = 1;
st(67,3,14) = -1;
st(7,4,14) = -1;
st(68,3,14) = 1;
st(68,3,15) = -1;
st(12,3,15) = -1;
st(69,3,15) = 1;
st(13,3,15) = 1;
st(69,3,16) = -1;
st(67,3,16) = 1;
st(8,4,16) = 1;
st(93,3,17) = -1;
st(13,3,17) = -1;
st(22,2,17) = -1;
st(94,3,17) = 1;
st(21,2,17) = 1;
st(27,2,17) = 2;
st(103,2,17) = 1;
st(94,3,18) = -1;
st(13,3,18) = -1;
st(22,2,18) = -1;
st(27,4,18) = -2;
st(103,2,18) = -1;
st(93,3,18) = 1;
st(12,3,18) = 1;
st(13,3,18) = 1;
st(21,3,18) = 1;
st(27,2,18) = 2;
st(95,3,19) = -1;
st(21,2,19) = -2;
st(96,3,19) = 1;
st(96,3,20) = -1;
st(28,4,20) = -1;
st(97,3,20) = 1;
st(22,2,20) = 2;
st(103,2,20) = 2;
st(97,3,21) = -1;
st(21,2,21) = -1;
st(27,4,21) = -2;
st(103,4,21) = -1;
st(98,3,21) = 1;
st(22,2,21) = 1;
st(27,2,21) = 1;
st(98,3,22) = -1;
st(21,2,22) = -1;
st(27,4,22) = -2;
st(99,3,22) = 1;
st(22,2,22) = 1;
st(27,2,22) = 1;
st(103,4,22) = 1;
st(99,3,23) = -1;
st(27,4,23) = -2;
st(103,2,23) = -2;
st(100,3,23) = 1;
st(29,4,23) = 1;
st(27,2,23) = 1;
st(100,3,24) = -1;
st(27,4,24) = -2;
st(95,3,24) = 1;
st(29,4,24) = 1;
st(27,2,24) = 1;
st(109,3,25) = -1;
st(104,2,25) = -1;
st(110,3,25) = 1;
st(110,3,26) = -1;
st(104,2,26) = -1;
st(111,3,26) = 1;
st(111,3,27) = -1;
st(110,3,27) = 1;
st(104,4,27) = 1;
st(110,3,28) = -1;
st(109,3,28) = 1;
st(104,4,28) = 1;
st(112,3,29) = -1;
st(104,4,29) = -1;
st(17,2,29) = -3;
st(113,3,29) = 1;
st(113,3,30) = -1;
st(112,4,30) = 1;
st(104,2,30) = 1;
st(17,3,30) = 3;

st(:) = 0;


%% Intra-Compartment Reaction Rate Constants
% (n_intra_reactions, 2)
% column 1 is forward rate, column 2 is back rate
ki(1,1) = 1.0e1;
ki(2,1) = 1.0e1;
ki(3,1) = 1.0e1;
ki(4,1) = 1.0e1;
ki(5,1) = 1.0e1;
ki(6,1) = 1.0e1;
ki(7,1) = 1.0e1;
ki(8,1) = 1.0e1;
ki(9,1) = 1.0e1;
ki(10,1) = 1.0e1;
ki(11,1) = 1.0e1;
ki(12,1) = 1.0e1;
ki(13,1) = 1.0e1;
ki(14,1) = 1.0e1;
ki(15,1) = 1.0e1;
ki(16,1) = 1.0e1;
ki(17,1) = 1.0e1;
ki(18,1) = 1.0e1;
ki(19,1) = 1.0e1;
ki(20,1) = 1.0e1;
ki(21,1) = 1.0e1;
ki(22,1) = 1.0e1;
ki(23,1) = 1.0e1;
ki(24,1) = 1.0e1;
ki(25,1) = 1.0e1;
ki(26,1) = 1.0e1;
ki(27,1) = 1.0e1;
ki(28,1) = 1.0e1;
ki(29,1) = 1.0e1;
ki(30,1) = 1.0e1;
ki(31,1) = 1.0e1;
ki(32,1) = 1.0e1;
ki(33,1) = 1.0e1;
ki(34,1) = 1.0e1;
ki(35,1) = 4.4693e2;
ki(36,1) = 4.4693e2;
ki(37,1) = 6.7615e6;
ki(38,1) = 1.15473e5;
ki(39,1) = 3.474e3;
ki(40,1) = 1.905286e6;
ki(41,1) = 1.7314583e7;
ki(42,1) = 3.266667e6;
ki(43,1) = 9.7e5;
ki(44,1) = 5.47857e5;


%% Trans-Compartment Reaction Rates,
% (n_trans_reactions, 2)
% column 1 is forward rate, column 2 is back rate
kt(1,1) = 1.22e5;
kt(2,1) = 1.0e11;
kt(3,1) = 1.0e11;
kt(4,1) = 1.22e5;
kt(5,1) = 1.0e11;
kt(6,1) = 6.889e5;
kt(7,1) = 1.0e11;
kt(8,1) = 1.0e11;
kt(9,1) = 3.7453;
kt(10,1) = 1.0e11;
kt(11,1) = 1.0e11;
kt(12,1) = 9.2083e6;
kt(13,1) = 1.0e11;
kt(14,1) = 5.567e5;
kt(15,1) = 1.0e11;
kt(16,1) = 1.0e11;
kt(17,1) = 1.0e11;
kt(18,1) = 1.03e6;
kt(19,1) = 1.0e11;
kt(20,1) = 6.0e9;
kt(21,1) = 1.0e11;
kt(22,1) = 1.0e11;
kt(23,1) = 1.0e11;
kt(24,1) = 1.0e11;
kt(25,1) = 1.0e11;
kt(26,1) = 1.0e11;
kt(27,1) = 1.0e5;
kt(28,1) = 1.0e11;
kt(29,1) = 1.0e11;
kt(30,1) = 1.0e5;
kt(:,2) = 1.0e-2;


%% concentration inside c(n_comp x n_species)
% this should actually be n_species, n_comp, but already writen
% as n_comp, n_species, so in order to save time, 
% we'll just take the transpose of it:)

% first set each concentration to a very small value, can not be zero.
c0(:) = 1.0e-10;

c0(18,4) = 2.3000052e-3;
c0(19,4) = 1.87e-4;
c0(27,4) = 5.20e-9;
c0(104,4) = 1.15e-3;
c0(16,2) = 9.6880398e-3;
c0(17,2) = 3.40570052e-2;
c0(18,2) = 1.87e-4;
c0(21,2) = 2.05e-3;
c0(22,2) = 2.73333333333333e-3;
c0(23,2) = 3.535e-3;
c0(24,2) = 1.7e-3;
c0(25,2) = 1.51e-4;
c0(27,2) = 3.98e-8;
c0(30,2) = 7.07e-3;
c0(31,2) = 5.1e-3;
c0(103,2) = 5.2e-9;
c0(104,2) = 1.15e-3;
c0(114,2) = 5.20417042793042e-18;
c0(12,3) = 1.9e-2;
c0(13,3) = 1.0e-3;
c0(67,3) = 8.9e-5;
c0(68,3) = 1.11095e-3;
c0(69,3) = 1.525e-4;
c0(79,3) = 1.0e-4;
c0(80,3) = 1.000052e-4;
c0(81,3) = 2.870052e-4;
c0(82,3) = 2.87e-4;
c0(83,3) = 3.17e-4;
c0(84,3) = 3.6935e-3;
c0(85,3) = 3.4875e-3;
c0(86,3) = 2.27175e-3;
c0(87,3) = 1.0e-4;
c0(88,3) = 7.96715788571429e-3;
c0(89,3) = 5.29e-5;
c0(90,3) = 1.0529e-3;
c0(91,3) = 2.00529e-2;
c0(92,3) = 2.0529e-3;
c0(93,3) = 1.11e-4;
c0(94,3) = 1.9111e-2;
c0(95,3) = 3.17e-4;
c0(96,3) = 1.62330237e-2;
c0(97,3) = 8.517e-3;
c0(98,3) = 8.517e-3;
c0(99,3) = 8.517e-3;
c0(100,3) = 8.517e-3;
c0(109,3) = 1.0e-4;
c0(110,3) = 6.25e-4;
c0(111,3) = 6.0e-4;
c0(112,3) = 1.0e-4;
c0(113,3) = 7.06140104e-3;
c0(1,4) = 3.9e-6;
c0(2,4) = 3.75e-4;
c0(3,4) = 3.75e-4;
c0(4,4) = 2.9e-5;
c0(5,4) = 2.02e-4;
c0(6,4) = 6.3e-4;
c0(7,4) = 1.065e-3;
c0(8,4) = 1.08e-4;
c0(9,4) = 4.91e-4;
c0(10,4) = 1.0e-5;
c0(11,4) = 1.8e-4;
c0(14,4) = 2.0e-3;
c0(15,4) = 1.0e-3;
c0(16,4) = 1.014060052e-1;
c0(17,4) = 1.860538e-1;
c0(18,4) = 1.87e-4;
c0(19,4) = 1.87e-4;
c0(20,4) = 2.14e-4;
c0(23,4) = 3.535e-3;
c0(24,4) = 1.7e-3;
c0(25,4) = 1.51e-4;
c0(26,4) = 4.863e-3;
c0(27,4) = 5.2e-9;
c0(28,4) = 1.0e-6;
c0(29,4) = 6.435e-1;
c0(30,4) = 3.535e-3;
c0(31,4) = 1.7e-3;
c0(32,4) = 1.07e-4;
c0(33,4) = 2.7e-5;
c0(34,4) = 2.7e-5;
c0(35,4) = 2.7e-5;
c0(36,4) = 4.28e-4;
c0(37,4) = 2.7e-5;
c0(38,4) = 2.7e-5;
c0(40,4) = 1.0e-4;
c0(41,4) = 7.17e-3;
c0(42,4) = 9e-3;
c0(43,4) = 1.725e-3;
c0(44,4) = 1.868e-3;
c0(45,4) = 1.1e-5;
c0(46,4) = 1.05e-5;
c0(47,4) = 1.245e-5;
c0(48,4) = 3.78666666666667e-4;
c0(49,4) = 1e-4;
c0(50,4) = 4.08333333333333e-4;
c0(51,4) = 4.08333333333333e-4;
c0(52,4) = 4.7e-5;
c0(53,4) = 2.047e-3;
c0(54,4) = 2.134e-3;
c0(55,4) = 1.047e-3;
c0(56,4) = 2.2e-5;
c0(57,4) = 2.322e-3;
c0(58,4) = 2.2e-5;
c0(59,4) = 1.076e-3;
c0(60,4) = 2.2e-5;
c0(61,4) = 2.2e-5;
c0(62,4) = 2.2e-5;
c0(63,4) = 1.08e-4;
c0(64,4) = 6.84e-4;
c0(65,4) = 2.077e-3;
c0(66,4) = 4.971e-3;
c0(70,4) = 5e-4;
c0(71,4) = 3.58e-4;
c0(72,4) = 7.41e-4;
c0(73,4) = 9.1e-5;
c0(74,4) = 2.091e-3;
c0(75,4) = 3.073e-3;
c0(76,4) = 1.1110052e-3;
c0(77,4) = 1.0e-3;
c0(78,4) = 1.0910052e-3;
c0(101,4) = 1.8e-4;
c0(102,4) = 2.087e-3;
c0(104,4) = 1.15e-3;
c0(105,4) = 1.161e-3;
c0(106,4) = 2.7e-5;
c0(107,4) = 1.161e-3;
c0(108,4) = 4.89e-3;
c0(114,4) = 1.202005498765e-2;
c0(115,4) = 4.668955e-3;

% voltage v(n_comp)
% initial values for compartment potentials. 
v0(:) = 0;


verify_stochiometry(si, z);

%% Run the simulation

% load the initial conditions from a file
% this reads a given filename, and populates the 
% variables 'c0' and 'v0' with the values stored 
% in the file. 
% notes: the mat file must contain the variables c0 and v0.
% load overwrites whatever the current value of c0 and v0
% with the values stored in the file. 
load('newlonger.mat');

% make the function that the integrator calls. 
fun = odefun(cap, a, l, h, z, o, si, ki, st, kt, r);

% pack the initial values into the state vector
state = karyote_pack(c0,v0);

yp0 = zeros(size(state));

% test the function, call it once with the starting state vector. 
%options = odeset('NonNegative', 1:length(state)-n_comp, ...
%        'MaxStep', 1e-5);
    
% test the function, call it once with the starting state vector. 
t0 = 0;
tf = 0.000001;

%options = odeset('NonNegative', 1:(length(state)-n_comp), ...
%                 'RelTol', 1e-15, ...
%                 'AbsTol', 1e-30, ...
%                 'InitialStep', 0.01*abs(t0-tf));
%[t,y] = ode15s(fun, [t0 tf], state, options);

options = odeset('NonNegative', 1:(n_species*n_comp), ...
                 'InitialStep', 0.01*abs(t0-tf)); 
[t,y] = ode15s(fun, [t0, tf], state, options);

c = y(:,1:(n_species*n_comp));
v = y(:,(n_species*n_comp)+1:end);

disp(c(1,:));
disp(c(end,:));
    
subplot(2,1,1);
plot(t,c);

subplot(2,1,2);
plot(t,v);