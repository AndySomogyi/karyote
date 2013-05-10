%% Test diffusion of different charged species
% Expected results is all species diffuse to zero, however
% there is not enough species in the compartment to bring the
% potential down as the capacitance is large. 

% set the global constants
global F R T
F= 9.64853399e4; R=8.3144621; T=300;

n_comp = 2;
n_species = 3;
n_reactions = 0;
n_trans_reactions = 0;

% valences, each species has same valence everywhere, n_species.
z = [-1, 2, 0];

% concentration, n_species, n_comp
% each column is a compartment, each row in the column is a species.
c0 = ones(n_species,n_comp) * 1e-15;
c0(1,2) = 2;
c0(2,2) = 1;
c0(3,2) = 3;

c0(:) = c0(:) * 1e-5;

c0 = neutralize_charge(c0,z);

disp('concentratoin after neutralization');
disp(c0);

% stoichiometry matrix n_species, n_comp, n_reactions
si = zeros(n_species, n_comp, n_reactions);
st = zeros(n_species, n_comp, n_trans_reactions);

% forward and backward rate constants, 
% n_reactions, 2.
ki = zeros(n_reactions,2);
kt = zeros(n_trans_reactions,2);

% initial compartment voltage
v = [0, .001];

cap = abs(eye(n_comp,n_comp) - 1) * 1e2;

a = abs(eye(n_comp,n_comp) - 1);
l = abs(eye(n_comp,n_comp) - 1);
o = ones(n_comp, 1);
o(1) = inf;

r = inf(n_comp,n_comp);

% membrane permeability
% test with membrane partially permeable to all species.
h = ones(n_species, n_comp, n_comp);

fun = odefun(cap, a, l, h, z, o, si, ki, st, kt, r);

state = karyote_pack(c0,v);

options = odeset('NonNegative', 1:(n_species*n_comp)); 
[t,y] = ode45(fun, [0 7], state, options);

c = y(:,1:(n_species*n_comp));
v = y(:,(n_species*n_comp)+1:end);

disp(c(1,:));
disp(c(end,:));
    
subplot(2,1,1);
plot(t,c);

subplot(2,1,2);
plot(t,v);