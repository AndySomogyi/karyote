function [dvdt] = dvdt(a, z, j, dvdt_inv, v, r)
% An eqn for each compartment, a is define by
% F \sum_{a'!=a} A_{a,a'} z \cdot J_{a,a'} = 
% \sum_{a'!=a} A_{a,a'}C_{a,a'} \frac{dV_{a,a'}{dt}
% The dv/dt terms can be solved by treating this as a linear system, 
% Ax=b, where x is the column vector of dv/dt and solving for it. 
%
% The convention is that the 1st compartment is grounded, so both
% v and dv/dt are zero, so it is not considered a state variable. 
%
% a: area of membrane separating compartment a from a'
% z: valence of species i
% j: intra compartment flux, n_species,n_comp,n_comp

% v: potential of each compartment, including the ground compartment, 
%    size: n_comp,1
% r: resistance between compartment, size: n_comp,n_comp.
%
% returns a row vector of dvdt
%
% convention is all capacitor currents are entering the compartment, 
% and all active process and resitive currents are leaving.
%
    global F 
    n = size(a,1);
    
    rhs = zeros(n,1);
    
    % net ionic current from active processes.
    zj = zeros(1,n);
    
    % resitive current leaving each compartment
    jr = zeros(size(rhs));
    
    % loop over compartments, 1 is ground comp.
    for i = 2:n
        for k = 1:n
            % ionic currents leaving the compartment
            zj(k) = z * (j(:,i,k) - j(:,k,i));
        end
        rhs(i) = F * a(i,:) * zj';
        
        % the resitive currents leaving the compartment
        jr(i) = sum((v(i)-v) ./ r(:,i));
    end
    
    rhs = rhs + jr;

    dvdt = [0; dvdt_inv * rhs(2:end)];
end
