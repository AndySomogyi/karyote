%% Test of capacitor charging with a constant current source
% Expected results is that the voltage initially rises
% in accordance with standard capacitor charging, but then drops as
% the compartment is filled and some current flows
% backwards across the membrane. 

% set the global constants
global F R T
F= 9.64853399e4; R=8.3144621; T=300;

n_comp = 2;
n_species = 2;
n_reactions = 0;
n_trans_reactions = 1;

% valences, each species has same valence everywhere, n_species.
z = [-1, 1];

% concentration, n_species, n_comp
% each column is a compartment, each row in the column is a species.
c0 = ones(n_species,n_comp) * 1e-1;
%c0(2,2) = 0.1010;

c0 = neutralize_charge(c0,z);

disp('concentration after neutralization');
disp(c0);

% intra stoichiometry matrix n_species, n_comp, n_reactions
si = zeros(n_species, n_comp, n_reactions);

% trans stoichiometry matrix n_species, n_comp, n_reactions
st = zeros(n_species, n_comp, n_trans_reactions);
st(2,1,1) =  -1;
st(2,2,1) = 1;

% forward and backward rate constants, 
% n_reactions, 2.
ki = zeros(n_reactions,2);

kt = zeros(n_trans_reactions,2);
kt(1,1) = 0.001;
kt(1,2) = 0;

% initial compartment voltage
v = [ 0    0];

cap = abs(eye(n_comp,n_comp) - 1) * 1e3;

a = abs(eye(n_comp,n_comp) - 1);
l = abs(eye(n_comp,n_comp) - 1);
o = ones(n_comp, 1);
o(2) = 0.0001;
o(1) = inf;

r = inf(n_comp,n_comp);

% membrane permeability
% test with membrane partially permeable to all species.
h = zeros(n_species, n_comp, n_comp);
h(1,:,:) = 0.01;
h(2,:,:) = 1;

fun = odefun(cap, a, l, h, z, o, si, ki, st, kt, r);

state = karyote_pack(c0,v);

t0 = 0;
tf = 0.025;
options = odeset('NonNegative', 1:(n_species*n_comp), ....
                 'MaxStep', 0.01*abs(t0-tf)); 
[t,y] = ode15s(fun, [t0 tf], state, options);

c = y(:,1:(n_species*n_comp));
v = y(:,(n_species*n_comp)+1:end);

disp(c(1,:));
disp(c(end,:));
    
subplot(2,1,1);
plot(t,c);

subplot(2,1,2);
plot(t,v);