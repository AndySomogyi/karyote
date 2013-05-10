% set the global constants
global F R T
F= 9.64853399e4; R=8.3144621; T=300;

% 
% 

n_comp = 1;
n_species = 10;
n_reactions = 1;
n_trans_reactions = 1;

% concentration, n_species, n_comp
% each column is a compartment, each row in the column is a species.
c0 = rand(n_species,n_comp);



  
% stoichiometry matrix n_species, n_comp, n_reactions
si = zeros(n_species, n_comp, n_reactions);
st = zeros(n_species, n_comp, n_trans_reactions);

% forward and backward rate constants, 
% n_reactions, 2.
ki = zeros(n_reactions,2);
kt = zeros(n_trans_reactions,2);


% valences, each species has same valence everywhere, n_species.
% z must be a row vector
z = rand(1,n_species)-.5;


% initial compartment voltage
v = [0, 1];

cap = abs(eye(n_comp,n_comp) - 1);

a = abs(eye(n_comp,n_comp) - 1);
l = abs(eye(n_comp,n_comp) - 1);
o = ones(n_comp, 1);
o(1) = inf;

r = inf(n_comp,n_comp);

% membrane permeability
% test with membrane partially permeable to all species.
h = ones(n_species, n_comp, n_comp) * 1;

%fun = odefun(cap, a, l, h, z, o, si, ki, st, kt, r);



c2 = neutralize_charge(c0, z);
c2 = neutralize_charge(c2, z);


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


