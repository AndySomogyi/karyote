% set the global constants
global F R T
F= 96485.3365; R=8.3144621; T=300;

% H2(g) + 2 ICl(g) -> I2(g) + 2 HCl(g)
% 1       1           3       4
n_comp = 2;
n_species = 1;
n_reactions = 0;

% concentration, n_species, n_comp
% each column is a compartment, each row in the column is a species.
c0 = [0;
      1]';
  
% stoichiometry matrix n_species, n_comp, n_reactions
s = zeros(n_species, n_comp, n_reactions);

% forward and backward rate constants, 
% n_reactions, 2.
k = [];

% valences, each species has same valence everywhere, n_species.
z = [-1];

% initial compartment voltage
v = [-.2 .2];

cap = abs(eye(n_comp,n_comp) - 1) * 2e-1;

a = abs(eye(n_comp,n_comp) - 1);
l = abs(eye(n_comp,n_comp) - 1);
o = ones(n_comp, 1);

% membrane permeability
h = zeros(n_species, n_comp, n_comp);
h(1,1,2) = 1;
h(1,2,1) = 1;
h(:) =1;
%h(3,1,2) = 1;
%h(3,2,1) = 1;

fun = odefun(cap, a, l, h, z, o, s, k);

state = karyote_pack(c0,v);

[t,y] = ode45(fun, [0 2], state);

c = y(:,1:(n_species*n_comp));
v = y(:,(n_species*n_comp)+1:end);


disp(c(1,:));
disp(c(end,:));
    

subplot(2,1,1);
plot(t,c);

subplot(2,1,2);
plot(t,v);


