%% Get mitochondira concentrations working.
% set the global constants
global F R T
F= 9.64853399e4; R=8.3144621; T=300;

% 
% 

n_comp = 4;
n_species = 1;
n_reactions = 0;
n_trans_reactions = 0;

% concentration, n_species, n_comp
% each column is a compartment, each row in the column is a species.
c0 = [0 1 0 0]  * 5.2e-6;

  
% stoichiometry matrix n_species, n_comp, n_reactions
si = zeros(n_species, n_comp, n_reactions);
st = zeros(n_species, n_comp, n_trans_reactions);

% forward and backward rate constants, 
% n_reactions, 2.
ki = zeros(n_reactions,2);
kt = zeros(n_trans_reactions,2);


% valences, each species has same valence everywhere, n_species.
z = [1];
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

% initial compartment voltage
v = [0, 1, 1, 1];

cap = abs(eye(n_comp,n_comp) - 1);

a = abs(eye(n_comp,n_comp) - 1);
l = abs(eye(n_comp,n_comp) - 1);
o = ones(n_comp, 1);
o(1) = inf;

r = inf(n_comp,n_comp);

% membrane permeability
% test with membrane partially permeable to all species.
h = ones(n_species, n_comp, n_comp) * 1;

fun = odefun(cap, a, l, h, z, o, si, ki, st, kt, r);

state = karyote_pack(c0,v);

[t,y] = ode15s(fun, [0 2], state);

c = y(:,1:(n_species*n_comp));
v = y(:,(n_species*n_comp)+1:end);


disp(c(1,:));
disp(c(end,:));
    
subplot(2,1,1);
plot(t,c);

subplot(2,1,2);
plot(t,v);


