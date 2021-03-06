% set the global constants
global F R T
F= 96485.3365; R=8.3144621; T=300;

% H2(g) + 2 ICl(g) -> I2(g) + 2 HCl(g)
% 1       1           3       4
n_comp = 2;
n_species = 4;
n_reactions = 1;

% concentration, n_species, n_comp
% each column is a compartment, each row in the column is a species.
c0 = [1 2 0.00001 0.000001;
      1 2 0.00001 0.000001]';
c0(:,1) = c0(:,1) * 2;
  
% stoichiometry matrix n_species, n_comp, n_reactions
s = zeros(n_species, n_comp, n_reactions);
s(:,1,1) = [-1 -2 1 2];
s(:,2,1) = [-1 -2 1 2];

% forward and backward rate constants, 
% n_reactions, 2.
k = [0.1, 0.00001];

% valences, each species has same valence everywhere, n_species.
z = [1 -1 1 -1];

% initial compartment voltage
v = [0, 1];

cap = abs(eye(n_comp,n_comp) - 1) * 0.1;

a = abs(eye(n_comp,n_comp) - 1);
l = abs(eye(n_comp,n_comp) - 1);
o = ones(n_comp, 1);

r = inf(n_comp,n_comp);

% membrane permeability
% test with membrane partially permeable to all species.
h = ones(n_species, n_comp, n_comp) * 0.000001;

fun = odefun(cap, a, l, h, z, o, s, k, r);

state = karyote_pack(c0,v);

% note, starting to be a stiff system, 
[t,y] = ode15s(fun, [0 10], state);

c = y(:,1:(n_species*n_comp));
v = y(:,(n_species*n_comp)+1:end);


disp(c(1,:));
disp(c(end,:));
    
subplot(2,1,1);
plot(t,c);

subplot(2,1,2);
plot(t,v);


