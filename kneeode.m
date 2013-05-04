function kneeode
%KNEEODE  The "knee problem" with Nonnegativity constraints.

% Problem parameter
epsilon = 1e-6;

y0 = 0.5;
xspan = [0, 2*pi];
 
% Solve without imposing constraints 
%options = [];
%[x1,y1] = ode15s(@odefcn,xspan,y0,options);
 
% Impose nonnegativity constraint
options = odeset('NonNegative',1);
[x2,y2] = ode15s(@odefcn,xspan,y0,options);
 
%figure
%plot(x1,y1,'b.-',x2,y2,'g-')
plot(x2,y2,'b.-')
%axis([0,2,-1,1]);
title('The "knee problem"');
legend('No constraints','nonnegativity')
xlabel('x');
ylabel('solution y')

   function yp = odefcn(x,y)
      %yp = ((1 - x)*y - y^2)/epsilon;
      disp(y);
      yp = -0.5;
   end  
end  % kneeode

% The derivative function is defined within nested function odefcn. 
% The value of epsilon used in odefcn is obtained from the outer function:

