cap = abs(eye(3)-1);

a = zeros(3);
a(1,2) = 1;
a(2,3) = 1;
a(1,3) = 1;
% make symmetric
a = a' + triu(a,1);

z = 1;
j = zeros(1,3,3);
j(1,2,1) = 1;
j(1,3,1) = 1;
global F;

F = 1;

C = a .* cap;
dvdt_inv = C(2:end,2:end);
for i = 1:size(dvdt_inv, 1)
    dvdt_inv(i,i) = -sum(C(i+1,:));
end
clear C;
disp(dvdt_inv);
dvdt_inv = inv(dvdt_inv);



dv = dvdt(a, z, j, dvdt_inv);

disp(dv);