function  test_active_flux()
% H2(g) + 2 ICl(g) -> I2(g) + 2 HCl(g)
% 1       1           3       4


    c0 = [1 2 0 0 ]';

    s = [-1; -2; 1; 2];

    k = [0.1, 0.0001];

    active_flux(c0,s,k)

    function dcdt = fun(~,c)
        dcdt = active_flux(c,s,k);
    end

    [t,y] = ode45(@fun, [0 500], c0);
    
    plot(t,y);

end

