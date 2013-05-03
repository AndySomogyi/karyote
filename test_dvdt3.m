function [] = test_dvdt3()
% A 3 compartment system, equiv of a two coupled RC circuits. 

    n_comp = 3;

    %cap = abs(eye(n_comp)-1) * 0.1;
    cap = zeros(n_comp);
    cap(1,2) = 0.1;
    cap(2,3) = 0.1;
    cap = cap + triu(cap)';

    a = zeros(n_comp);
    a(1,2) = 1;
    a(2,3) = 1;
    a(1,3) = 1;
    % make symmetric
    a = a' + triu(a,1);

    z = 1;
    j = zeros(1,n_comp, n_comp);
    j(1,1,2) = 1;
    j(1,3,1) = 1;

    % resistance:
    r = inf(n_comp);
    r(1,2) = 1;
    r(2,3) = 1;
    r(1,3) = 1;
    r = triu(r) + triu(r,1)';
    global F;

    F = 1;

    % convention is all capacitor currents are entering the compartment, 
    % and all active process and resitive currents are leaving. 
    C = a .* cap;
    dvdt_inv = C(2:end,2:end);
    for i = 1:size(dvdt_inv, 1)
        dvdt_inv(i,i) = -sum(C(i+1,:));
    end
    clear C;
    disp(dvdt_inv);
    dvdt_inv = inv(dvdt_inv);
    assert(all(all(isfinite(dvdt_inv))), 'capacitance matrix is not invertable');

    v0 = [0 -1 5]';



    dv = dvdt(a, z, j, dvdt_inv, v0, r);
    
    disp(dv);


    function [state] = fun(t,state)
        state = dvdt(a,z,j,dvdt_inv, state, r);
    end

    [t,y] = ode45(@fun,[0 1], v0);
    
    plot(t,y);

end
