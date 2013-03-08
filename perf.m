function [t,y] = perf(a)
[f,init] = test(5,10,15);
[t,y] = ode45(f,[0 0.1], init);
end
