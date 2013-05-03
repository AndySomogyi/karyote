function [t,y] = test_dvdt1()

    cap = abs(eye(3)-1) * 1;

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

    % resistance:
    r = ones(3,3) * 1000;
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

    v0 = [0 1 2]';



    dv = dvdt(a, z, j, dvdt_inv, v0, r);
    
    disp(dv);


    function [state] = fun(t,state)
        state = dvdt(a,z,j,dvdt_inv, state, r);
    end

    [t,y] = ode45(@fun,[0 1], v0);
    
    plot(t,y);

    

end
