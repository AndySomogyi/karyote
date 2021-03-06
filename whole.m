% set the global constants
global F R T
F= 96485.3365; R=8.3144621; T=300;

% define the size of the system
% the semicolon at the end is not required, all it does is suppress 
% diplaying of the value..
% compartment 1 is ground and 2 is matrix
n_comp = 3;
n_species = 85;
n_intra_reactions = 47;
n_trans_reactions = 9;

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
a(1, 3) = 1.866e-4;
a(2, 3) = 1.866e-4;

% make symmetric
a = a' + triu(a,1);

% thickness of membrane in nm
% default is infinitly thick membrane
l(:) = Inf;

l(3,1) = 4;
l(2,3) = 2;

l(1,3) = 4;
l(3,2) = 2;

%% permeability, (n_species, n_comp, n_comp) %%
% diffuse from universe to into outer comp (1)
h(72,1,3) = 1;
h(72,3,1) = 1;
h(85,2,3) = 1;
h(85,3,2) = 1;

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
z(21) = -4;
z(22) = -3;
z(23) = -1;
z(24) =  1;
z(26) = -2;
z(27) = -1;
z(28) = -2;
z(37) = -2;
z(38) = -2;
z(39) = -3;
z(41) = -3;
z(42) = -3;
z(44) =  1;
z(46) =  4;
z(50) = -2;
z(55) = -2;
z(56) = -3;
z(57) = -1;
z(59) = -2;
z(60) = -2;
z(62) = -2;
z(63) = -2;
z(65) =  1;
z(66) = -1;
z(67) = -1;
z(68) =  2;
z(69) =  1;
z(70) = -1;
z(71) = -1;
z(72) =  2;
z(73) =  2;
z(75) =  2;
z(76) = -1;
z(78) = -1;
z(79) = -2;
z(80) = -4;
z(81) = -4;
z(82) = -2;
z(83) =  2;
z(84) =  3;


%% volume of each compartment
o(:) = [Inf, 1.17e-6 1.6e-2];

%% stoichometry of intra-compartment reactions, indexing: nspecies x ncomp x nreactions
si(74,2,1) = -1;
si(72,2,1) = -1;
si(75,2,1) =  1;
si(75,2,2) = -1;
si(76,2,2) = -1;
si(74,2,2) = 1;
si(23,2,2) = 1;
si(29,2,2) = 1;
si(72,2,2) = 1;
si(18,2,3) = -1;
si(29,2,3) = -1;
si(24,2,3) = -1;
si(20,2,3) = 1;
si(30,2,3) = 1;
si(30,2,4) = -1;
si(31,2,4) = -1;
si(29,2,4) =  1;
si(32,2,4) =  1;
si(32,2,5) = -1;
si(11,2,5) = -1;
si(33,2,5) =  1;
si(1,2,5) =  1;
si(33,2,6) = -1;
si(35,2,6) = -1;
si(31,2,6) =  1;
si(34,2,6) =  1;
si(34,2,7) = -1;
si(14,2,7) = -1;
si(35,2,7) =  1;
si(15,2,7) =  1;
si(24,2,7) =  1;
si(36,2,8) = -1;
si(10,2,8) = -1;
si(37,2,8) =  1;
si(37,2,9) = -1;
si(1,2,9) = -1;
si(38,2,9) = 1;
si(38,2,10) = -1;
si(39,2,10) = 1;
si(70,2,10) = 1;
si(24,2,10) = 2;
si(39,2,11) = -1;
si(36,2,11) = 1;
si(2,2,11) = 1;
si(40,2,12) = -1;
si(2,2,12) = -1;
si(41,2,12) = 1;
si(41,2,13) = -1;
si(42,2,13) = 1;
si(25,2,13) = 1;
si(42,2,14) = -1;
si(25,2,14) = -1;
si(40,2,14) = 1;
si(4,2,14) = 1;
si(43,2,15) = -1;
si(14,2,15) = -1;
si(44,2,15) =  1;
si(44,2,16) = -1;
si(71,2,16) = -1;
si(45,2,16) = 1;
si(45,2,17) = -1;
si(72,2,17) = -2;
si(46,2,17) = 1;
si(5,2,17) = 1;
si(28,2,17) = 1;
si(24,2,17) = 2;
si(68,2,17) = 1;
si(46,2,18) = -1;
si(43,2,18) = 1;
si(15,2,18) = 1;
si(72,2,18) = 2;
si(47,2,19) = -1;
si(72,2,19) = -1;
si(73,2,19) = 1;
si(5,2,20) = -1;
si(73,2,20) = -1;
si(48,2,20) = 1;
si(20,2,20) = 1;
si(48,2,21) = -1;
si(49,2,21) = -1;
si(73,2,21) = 1;
si(50,2,21) = 1;
si(50,2,22) = -1;
si(11,2,22) = -1;
si(51,2,22) = 1;
si(6,2,22) = 1;
si(51,2,23) = -1;
si(53,2,23) = -1;
si(49,2,23) = 1;
si(52,2,23) = 1;
si(52,2,24) = -1;
si(14,2,24) = -1;
si(53,2,24) = 1;
si(15,2,24) = 1;
si(24,2,24) = 1;
si(54,2,25) = -1;
si(6,2,25) = -1;
si(55,2,25) = 1;
si(55,2,26) = -1;
si(23,2,26) = -1;
si(56,2,26) = 1;
si(11,2,26) = 1;
si(56,2,27) = -1;
si(57,2,27) = 1;
si(7,2,27) = 1;
si(57,2,28) = -1;
si(22,2,28) = -1;
si(54,2,28) = 1;
si(21,2,28) = 1;
si(58,2,29) = -1;
si(7,2,29) = -1;
si(59,2,29) = 1;
si(59,2,30) = -1;
si(12,2,30) = -1;
si(60,2,30) = 1;
si(13,2,30) = 1;
si(60,2,31) = -1;
si(58,2,31) = 1;
si(8,2,31) = 1;
si(61,2,32) = -1;
si(8,2,32) = -1;
si(62,2,32) = 1;
si(62,2,33) = -1;
si(25,2,33) = -1;
si(63,2,33) = 1;
si(63,2,34) = -1;
si(61,2,34) = 1;
si(9,2,34) = 1;
si(64,2,35) = -1;
si(14,2,35) = -1;
si(65,2,35) = 1;
si(65,2,36) = -1;
si(9,2,36) = -1;
si(66,2,36) = 1;
si(66,2,37) = -1;
si(67,2,37) = 1;
si(67,2,38) = -1;
si(69,2,38) = 1;
si(10,2,38) = 1;
si(69,2,39) = -1;
si(64,2,39) = 1;
si(15,2,39) = 1;
si(24,2,39) = 1;
si(19,2,40) = -1;
si(18,2,40) = 1;
si(24,2,40) = 1;
si(21,2,41) = -1;
si(68,2,41) = -1;
si(26,2,41) = 1;
si(22,2,42) = -1;
si(68,2,42) = -1;
si(27,2,42) = 1;
si(77,2,43) = -1;
si(27,2,43) = -1;
si(78,2,43) =  1;
si(78,2,44) = -1;
si(22,2,44) = -1;
si(80,2,44) =  1;
si(80,2,45) = -1;
si(81,2,45) = 1;
si(81,2,46) = -1;
si(79,2,46) =  1;
si(82,2,46) =  1;
si(79,2,47) = -1;
si(77,2,47) =  1;
si(26,2,47) =  1;

%% stoichometry of transmembrane reactions, indexing: nspecies x ncomp x nreactions 
st(18,1,1) = -1;
st(24,1,1) = -1;
st(18,2,1) =  1;
st(24,2,1) =  1;
st(72,1,2) = -2;
st(72,2,2) =  2;
st(27,2,3) = -1;
st(23,2,3) = -1;
st(24,3,3) = -4;
st(26,2,3) =  1;
st(24,2,3) =  4;
st(22,3,4) = -1;
st(21,2,4) = -1;
st(22,2,4) =  1;
st(21,3,4) =  1;
st(15,2,5) = -1;
st(12,2,5) = -1;
st(24,2,5) = -5;
st(13,2,5) =  1;
st(24,3,5) =  4;
st(14,2,5) =  1;
st(13,2,6) = -1;
st(84,3,6) = -2;
st(24,2,6) = -2;
st(12,2,6) =  1;
st(83,3,6) =  1;
st(24,3,6) =  4;
st(83,3,7) = -4;
st(24,2,7) = -8;
st(85,2,7) = -1;
st(84,3,7) =  4;
st(25,2,7) =  2;
st(24,3,7) =  4;
st(72,2,8) = -1;
st(17,3,8) = -3;
st(72,3,8) =  1;
st(17,2,8) =  3;
st(72,2,9) = -2;
st(72,3,9) =  2;


%% Intra-Compartment Reaction Rate Constants
% (n_intra_reactions, 2)
% column 1 is forward rate, column 2 is back rate
ki(1,1) = 1.0e9;
ki(2,1) = 2.28e7;
ki(3,1) = 1.84e6;
ki(4,1) = 1.0e9;
ki(5,1) = 1.0e9;
ki(6,1) = 1.0e9;
ki(7,1) = 1.0e9;
ki(8,1) = 1.0e9;
ki(9,1) = 1.0e9;
ki(10,1) = 3.12e6;
ki(11,1) = 1.0e9;
ki(12,1) = 1.0e9;
ki(13,1) = 9.43e6;
ki(14,1) = 1.0e9;
ki(15,1) = 1.0e9;
ki(16,1) = 1.0e9;
ki(17,1) = 1.36e4;
ki(18,1) = 1.0e9;
ki(19,1) = 1.0e9;
ki(20,1) = 4.9e5;
ki(21,1) = 1.0e9;
ki(22,1) = 1.0e9;
ki(23,1) = 1.0e9;
ki(24,1) = 1.0e9;
ki(25,1) = 1.0e9;
ki(26,1) = 7.86e5;
ki(27,1) = 1.0e9;
ki(28,1) = 1.0e9;
ki(29,1) = 2.83e5;
ki(30,1) = 1.0e9;
ki(31,1) = 1.0e9;
ki(32,1) = 1.0e9;
ki(33,1) = 8.62e4;
ki(34,1) = 1.0e9;
ki(35,1) = 1.0e9;
ki(36,1) = 1.0e9;
ki(37,1) = 1.73e8;
ki(38,1) = 1.0e9;
ki(39,1) = 1.0e9;
ki(40,1) = 1.0e9;
ki(41,1) = 1.0e9;
ki(42,1) = 1.0e9;
ki(43,1) = 1.0e9;
ki(44,1) = 1.0e9;
ki(45,1) = 6.7615e6;
ki(46,1) = 1.0e9;
ki(47,1) = 1.0e9;

%% Trans-Compartment Reaction Rates,
% (n_trans_reactions, 2)
% column 1 is forward rate, column 2 is back rate
kt(1,1) = 9.23e3;
kt(2,1) = 2e1;
kt(3,1) = 4.76e3;
kt(4,1) = 6.76;
kt(5,1) = 7.34e5;
kt(6,1) = 1.03e6;
kt(7,1) = 4e7;
kt(8,1) = 1.0e1;
kt(9,1) = 1e3


%% concentration inside c(n_comp x n_species)
% this should actually be n_species, n_comp, but already writen
% as n_comp, n_species, so in order to save time, 
% we'll just take the transpose of it:)

% first set each concentration to a very small value, can not be zero.
% 1: ground, 2: matrix, 3: IMS
c0(18,1) = 2.3000052e-3;
c0(19,1) = 1.87e-4;
c0(24,1) = 5.20e-9;
c0(72,1) = 1.15e-3;
c0(1,2) = 3.9e-6;
c0(2,2) = 3.75e-4;
c0(3,2) = 3.75e-4;
c0(4,2) = 2.9e-5;
c0(5,2) = 2.02e-4;
c0(6,2) = 6.3e-4;
c0(7,2) = 1.065e-3;
c0(8,2) = 1.08e-4;
c0(9,2) = 4.91e-4;
c0(10,2) = 1.0e-5;
c0(11,2) = 1.8e-4;
c0(12,2) = 1.9e-2;
c0(13,2) = 1.0e-3;
c0(14,2) = 2.0e-3;
c0(15,2) = 1.0e-3;
c0(16,2) = 1.014060052e-1;
c0(17,2) = 1.860538e-1;
c0(18,2) = 1.87e-4;
c0(19,2) = 1.87e-4;
c0(20,2) = 2.14e-4;
c0(21,2) = 3.535e-3;
c0(22,2) = 1.7e-3;
c0(23,2) = 4.863e-3;
c0(24,2) = 5.2e-9;
c0(25,2) = 6.435e-1;
c0(26,2) = 3.535e-3;
c0(27,2) = 1.7e-3;
c0(28,2) = 1.07e-4;
c0(29,2) = 2.7e-5;
c0(30,2) = 2.7e-5;
c0(31,2) = 2.7e-5;
c0(32,2) = 4.28e-4;
c0(33,2) = 2.7e-5;
c0(34,2) = 2.7e-5;
c0(35,2) = 2.7e-5;
c0(36,2) = 1.1e-5;
c0(37,2) = 1.05e-5;
c0(38,2) = 1.245e-5;
c0(39,2) = 3.78666666666667e-4;
c0(40,2) = 1e-4;
c0(41,2) = 4.08333333333333e-4;
c0(42,2) = 4.08333333333333e-4;
c0(43,2) = 4.7e-5;
c0(44,2) = 2.047e-3;
c0(45,2) = 2.134e-3;
c0(46,2) = 1.047e-3;
c0(47,2) = 2.2e-5;
c0(48,2) = 2.322e-3;
c0(49,2) = 2.2e-5;
c0(50,2) = 1.076e-3;
c0(51,2) = 2.2e-5;
c0(52,2) = 2.2e-5;
c0(53,2) = 2.2e-5;
c0(54,2) = 1.08e-4;
c0(55,2) = 6.84e-4;
c0(56,2) = 2.077e-3;
c0(57,2) = 4.971e-3;
c0(58,2) = 8.9e-5;
c0(59,2) = 1.11095e-3;
c0(60,2) = 1.525e-4;
c0(61,2) = 5e-4;
c0(62,2) = 3.58e-4;
c0(63,2) = 7.41e-4;
c0(64,2) = 9.1e-5;
c0(65,2) = 2.091e-3;
c0(66,2) = 3.073e-3;
c0(67,2) = 1.1110052e-3;
c0(68,2) = 1.0e-3;
c0(69,2) = 1.0910052e-3;
c0(70,2) = 1.8e-4;
c0(71,2) = 2.087e-3;
c0(72,2) = 1.15e-10;
c0(73,2) = 1.161e-3;
c0(74,2) = 2.7e-5;
c0(75,2) = 1.161e-3;
c0(76,2) = 4.89e-3;
c0(77,2) = 1.0e-4;
c0(78,2) = 7.17e-3;
c0(79,2) = 9e-3;
c0(80,2) = 1.725e-3;
c0(81,2) = 1.868e-3;
c0(82,2) = 1.51e-4;
c0(16,3) = 9.6880398e-3;
c0(17,3) = 3.40570052e-2;
c0(18,3) = 1.87e-4;
c0(83,3) = 2.05e-3;
c0(84,3) = 2.73333333333333e-3;
c0(85,3) = 1.0e-6;
c0(21,3) = 3.535e-3;
c0(22,3) = 1.7e-3;
c0(82,3) = 1.51e-4;
c0(24,3) = 3.98e-8;
c0(26,3) = 7.07e-3;
c0(27,3) = 5.1e-3;
c0(72,3) = 1.15e-10;

c0 = neutralize_charge(c0, z);

% voltage v(n_comp)
% initial values for compartment potentials. 
v0(:) = 0;
v0(2) = .2795;
v0(3) = .2;

verify_stochiometry(si, z);

%% Run the simulation

% load the initial conditions from a file
% this reads a given filename, and populates the 
% variables 'c0' and 'v0' with the values stored 
% in the file. 
% notes: the mat file must contain the variables c0 and v0.
% load overwrites whatever the current value of c0 and v0
% with the values stored in the file. 

% make the function that the integrator calls. 
fun = odefun(cap, a, l, h, z, o, si, ki, st, kt, r);

% load initial conditions
load('whole.mat')


% pack the initial values into the state vector
state = karyote_pack(c0,v0);

% test the function, call it once with the starting state vector. 
%options = odeset('NonNegative', 1:length(state)-n_comp, ...
%        'MaxStep', 1e-5);
    
% test the function, call it once with the starting state vector. 
t0 = 0;
tf = 2e-8;

%options = odeset('NonNegative', 1:(length(state)-n_comp), ...
%                 'RelTol', 1e-15, ...
%                 'AbsTol', 1e-30, ...
%                 'InitialStep', 0.01*abs(t0-tf));
%[t,y] = ode15s(fun, [t0 tf], state, options);

options = odeset('NonNegative', 1:(n_species*n_comp), ...
                  'AbsTol', 1e-15, ...
                 'InitialStep', 0.00001*abs(t0-tf));
[t,y] = ode15s(fun, [t0, tf], state, options);

[c,v] = karyote_unpack(y, n_comp, n_species);

% c is now a [n_time, n_comp, n_species] matrix. 

% to display a single species, say 84 for all compartments, pick 
% that one out via:
squeeze(c(end,72,2))
% the squeeze() function is required because because MATALB is a little
% dumb when it comes to multi dim matricies: here we pick out the last 
% time value (end), all the compartments (:) and species 84 (84), so this
% results in a [1,3,1] matrix, MATALB is evidently not smart enough to
% recognize that this is a just a length 3 vector, so we have to tell
% is to squeeze it down to a length 3 vector via squeeze. 

% same thing with all the plots below:

disp('final values: ');
disp(squeeze(c(end,:,:)));

subplot(4,1,1);
plot(t,squeeze(c(:,:,1)));

subplot(4,1,2);
plot(t,squeeze(c(:,:,2)));

subplot(4,1,3);
plot(t,squeeze(c(:,:,3)));

subplot(4,1,4);
plot(t,v);