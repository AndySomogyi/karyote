% set the global constants
global F R T
F= 96485.3365; R=8.3144621; T=300;

% define the size of the system
% the semicolon at the end is not required, all it does is suppress 
% diplaying of the value..
n_comp = 2;
n_species = 115;
n_reactions = 34;

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

%    s: = stoichiometric coefficient of species i for reaction k in
%        compartment a. dims: nspecies, n_species x ncomp x n_reactions
s = zeros(n_species, n_comp, n_reactions);

%    k: rate (n_reaction, rate direction), 1 = forward, 2 = reverse
k = zeros(n_reactions, 2);

%    concentration inside c(n_comp x n_specie)
c = zeros(n_species, n_comp);

%    voltage
v = zeros(1,n_comp);

% "atom" count per species, not actually the real atom count, just
% how many things have to be conserved, for example, Adenesine is a 1
% because is always remains intact, but combines with other things.
ac = zeros(n_species);


% now we can start assigning actual values


% capacitance is easy, they are all the same values, 
cap(:) = 2e-04;

% area between compartments, 
a(1, 2) = 1.866e-4;

% make symmetric
a = a' + triu(a,1);

% thickness of membrane in nm
% default is infinitly thick membrane
l(:) = Inf;

% make symmetric
l = l' + triu(l,1);

% permeability, (n_comp, n_comp, n_species)

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
o(:) = [Inf 1.6e-04];

% stoichometry, indexing: nspecies x ncomp x nreactions
% 1-34 are fast rxns, 35-44 are slow rxns (fast and slow are intracomparmental) and 45-74 are
% transmembrane rxns


% make the reactions for what what previously comp1, 
comp1 = zeros(n_species, n_reactions);
comp1(106,1) = -1;
comp1(104,1) = -1;
comp1(107,1) =  1;
comp1(34,2) = -1;
comp1(35,2) = -1;
comp1(33,2) =  1;
comp1(36,2) =  1;
comp1(36,3) = -1;
comp1(11,3) = -1;
comp1(37,3) =  1;
comp1(1,3) =  1;
comp1(37,4) = -1;
comp1(39,4) = -1;
comp1(35,4) =  1;
comp1(38,4) =  1;
comp1(38,5) = -1;
comp1(14,5) = -1;
comp1(39,5) =  1;
comp1(15,5) =  1;
comp1(27,5) =  1;
comp1(40,6) = -1;
comp1(31,6) = -1;
comp1(41,6) =  1;
comp1(41,7) = -1;
comp1(24,7) = -1;
comp1(43,7) =  1;
comp1(44,8) = -1;
comp1(42,8) =  1;
comp1(25,8) =  1;
comp1(42,9) = -1;
comp1(40,9) =  1;
comp1(30,9) =  1;
comp1(45,10) = -1;
comp1(10,10) = -1;
comp1(46,10) =  1;
comp1(46,11) = -1;
comp1(1,11) = -1;
comp1(47,11) = 1;
comp1(48,12) = -1;
comp1(45,12) = 1;
comp1(2,12) = 1;
comp1(49,13) = -1;
comp1(2,13) = -1;
comp1(50,13) = 1;
comp1(51,14) = -1;
comp1(29,14) = -1;
comp1(49,14) = 1;
comp1(4,14) = 1;
comp1(52,15) = -1;
comp1(14,15) = -1;
comp1(53,15) =  1;
comp1(53,16) = -1;
comp1(102,16) = -1;
comp1(54,16) = 1;
comp1(55,17) = -1;
comp1(52,17) = 1;
comp1(15,17) = 1;
comp1(104,17) = 2;
comp1(56,18) = -1;
comp1(104,18) = -1;
comp1(105,18) = 1;
comp1(57,19) = -1;
comp1(58,19) = -1;
comp1(105,19) = 1;
comp1(59,19) = 1;
comp1(59,20) = -1;
comp1(11,20) = -1;
comp1(60,20) = 1;
comp1(6,20) = 1;
comp1(60,21) = -1;
comp1(62,21) = -1;
comp1(58,21) = 1;
comp1(61,21) = 1;
comp1(61,22) = -1;
comp1(14,22) = -1;
comp1(62,22) = 1;
comp1(15,22) = 1;
comp1(27,22) = 1;
comp1(63,23) = -1;
comp1(6,23) = -1;
comp1(64,23) = 1;
comp1(65,24) = -1;
comp1(66,24) = 1;
comp1(7,24) = 1;
comp1(66,25) = -1;
comp1(24,25) = -1;
comp1(63,25) = 1;
comp1(23,25) = 1;
comp1(70,26) = -1;
comp1(8,26) = -1;
comp1(71,26) = 1;
comp1(72,27) = -1;
comp1(70,27) = 1;
comp1(9,27) = 1;
comp1(73,28) = -1;
comp1(14,28) = -1;
comp1(74,28) = 1;
comp1(74,29) = -1;
comp1(9,29) = -1;
comp1(75,29) = 1;
comp1(76,30) = -1;
comp1(78,30) = 1;
comp1(10,30) = 1;
comp1(78,31) = -1;
comp1(73,31) = 1;
comp1(15,31) = 1;
comp1(27,31) = 1;

comp1(19,32) = -1;
comp1(18,32) = 1;
comp1(27,32) = 1;

comp1(23,33) = -1;
comp1(77,33) = -1;
comp1(30,33) = 1;
comp1(24,34) = -1;
comp1(77,34) = -1;
comp1(31,34) = 1;

% the world compartment, 
% no reactions go on in the world. 
world = zeros(n_species, n_reactions);

% populate the full stoichiometry matrix
s(:,1,:) = world;
s(:,2,:) = comp1;

%     k = (n_reactions, 2);

% the reverse rate, very small indicating irreversible reactions
k(:,2) = 1.0e-2;

% forward rate.
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


% concentration inside c(n_comp x n_specie)
% this should actually be n_species, n_comp, but already writen
% as n_comp, n_species, so in order to save time, 
% we'll just take the transpose of it:)

% first set each concentration to a very small value, can not be zero.
c(:) =  1.0e-10;

% comp1 
c(1,2) = 3.9e-6;
c(2,2) = 3.75e-4;
c(3,2) = 3.75e-4;
c(4,2) = 2.9e-5;
c(5,2) = 2.02e-4;
c(6,2) = 6.3e-4;
c(7,2) = 1.065e-3;
c(8,2) = 1.08e-4;
c(9,2) = 4.91e-4;
c(10,2) = 1.0e-5;
c(11,2) = 1.8e-4;
c(14,2) = 2.0e-3;
c(15,2) = 1.0e-3;
c(16,2) = 1.014060052e-1;
c(17,2) = 1.860538e-1;
c(18,2) = 1.87e-4;
c(19,2) = 1.87e-4;
c(20,2) = 2.14e-4;
c(23,2) = 3.535e-3;
c(24,2) = 1.7e-3;
c(25,2) = 1.51e-4;
c(26,2) = 4.863e-3;
c(27,2) = 5.2e-9;
c(28,2) = 1.0e-6;
c(29,2) = 6.435e-1;
c(30,2) = 3.535e-3;
c(31,2) = 1.7e-3;
c(32,2) = 1.07e-4;
c(33,2) = 2.7e-5;
c(34,2) = 2.7e-5;
c(35,2) = 2.7e-5;
c(36,2) = 4.28e-4;
c(37,2) = 2.7e-5;
c(38,2) = 2.7e-5;
c(40,2) = 1.0e-4;
c(41,2) = 7.17e-3;
c(42,2) = 9e-3;
c(43,2) = 1.725e-3;
c(44,2) = 1.868e-3;
c(45,2) = 1.1e-5;
c(46,2) = 1.05e-5;
c(47,2) = 1.245e-5;
c(48,2) = 3.78666666666667e-4;
c(49,2) = 1e-4;
c(50,2) = 4.08333333333333e-4;
c(51,2) = 4.08333333333333e-4;
c(52,2) = 4.7e-5;
c(53,2) = 2.047e-3;
c(54,2) = 2.134e-3;
c(55,2) = 1.047e-3;
c(56,2) = 2.2e-5;
c(57,2) = 2.322e-3;
c(58,2) = 2.2e-5;
c(59,2) = 1.076e-3;
c(60,2) = 2.2e-5;
c(61,2) = 2.2e-5;
c(62,2) = 2.2e-5;
c(63,2) = 1.08e-4;
c(64,2) = 6.84e-4;
c(65,2) = 2.077e-3;
c(66,2) = 4.971e-3;
c(70,2) = 5e-4;
c(71,2) = 3.58e-4;
c(72,2) = 7.41e-4;
c(73,2) = 9.1e-5;
c(74,2) = 2.091e-3;
c(75,2) = 3.073e-3;
c(76,2) = 1.1110052e-3;
c(77,2) = 1.0e-3;
c(78,2) = 1.0910052e-3;
c(101,2) = 1.8e-4;
c(102,2) = 2.087e-3;
c(104,2) = 1.15e-3;
c(105,2) = 1.161e-3;
c(106,2) = 2.7e-5;
c(107,2) = 1.161e-3;
c(108,2) = 4.89e-3;
c(114,2) = 1.202005498765e-2;
c(115,2) = 4.668955e-3;

% world concentration
c(18,1) = 2.3000052e-3;
c(19,1) = 1.87e-4;
c(27,1) = 5.20e-9;
c(104,1) = 1.15e-3;


% voltage v(n_comp)
% initial values for compartment potentials. 
v(:) = 0;

verify_stochiometry(s, z);

% make the function that the integrator calls. 
fun = odefun(cap, a, l, h, z, o, s, k);

% pack the initial values into the state vector
state = karyote_pack(c,v);

% test the function, call it once with the starting state vector. 
t0 = 0;
tf = 3.5;
options = odeset('NonNegative', 1:(length(state)-n_comp), ...
                 'RelTol', 1e-3, ...
                 'AbsTol', 1e-20, ...
                 'InitialStep', 0.001*abs(t0-tf));
[t,y] = ode45(fun, [t0 tf], state, options);

c = y(:,n_species+1:2*n_species);

disp(c(1,:));
disp(c(end,:));

v = y(:,end-3:end);

plot(t,c);









