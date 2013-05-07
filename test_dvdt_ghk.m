function [] = test_dvdt_ghk()
% A single compartment system, equiv of a single RC circuit. 

    n_comp = 2;

    cap = abs(eye(n_comp)-1) * 0.1;

    a = zeros(n_comp);
    a(1,2) = 1;
    % make symmetric
    a = a' + triu(a,1);

    z = 1;
    j = zeros(1,n_comp, n_comp);
    j(1,1,2) = 1;

    % resistance:
    r = ones(n_comp) * 1;
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

    v0 = [0 -1]';



    dv = dvdt(a, z, j, dvdt_inv, v0, r);
    
    disp(dv);


    function [state] = fun(t,state)
        state = dvdt(a,z,j,dvdt_inv, state, r);
    end

    %[t,y] = ode45(@fun,[0 1], v0);
    
    %plot(t,y);
    
    z = 1;
    v = [0 1];
    c = [0 1];
    j = ghk_flux(z,v(1) - v(2), c(1), c(2));
    disp(j);

end

function [j] = ghk_flux(z, v, c0, c1)   
    R = 1;
    T = 1;
    F = 1;
    
    if v ~= 0
        ev = exp((-v.* z) * F / R / T); 



        a = (((z * F).^2) .* v) / R * T;
        b = (c0 - c1 .* ev) ./ (1 - ev);
        j = -a .* b;
    else
    end
end
