function [dvdt] = dvdt(a, z, j, dvdt_inv, v, r)
% calculate dvdt for each compartment. 
%
% The convention is that the 1st compartment is grounded, so both
% v and dv/dt are zero, so it is not considered a state variable. 
%
% a: area of membrane separating compartment a from a', 
% size: n_comp,n_comp
% z: valence of species i, size: n_comp,1
% j: intra compartment flux, size: n_species,n_comp,n_comp
% dvdt_inv: inverse of capacitance matrix, size: (n_comp-1),(n_comp-1).
% v: potential of each compartment, including the ground compartment, 
%    size: n_comp,1
% r: resistance between compartment, size: n_comp,n_comp.
%
% returns a row vector of dvdt
%
% convention is all capacitor currents are entering the compartment, 
% and all active process and resitive currents are leaving.

    global F 
    n = size(a,1);
    
    rhs = zeros(n,1);
    
    % net ionic current from active processes.
    zj = zeros(1,n);
    
    % resitive current leaving each compartment
    jr = zeros(size(rhs));
    
    % loop over compartments, 1 is ground comp.
    % the 1st compartment is does not change, so we do not calculate
    % deriv for it. If we did, it would be a linear combination 
    % of the others and lead to a non-invertable syste. 
    % cleaner to start index at 2, and set a zero
    % at the last line rather than add 1 to get voltages and 
    % currents. 
    for i = 2:n
        for k = 1:n
            % ionic currents leaving the compartment
            zj(k) = z * (j(:,i,k) - j(:,k,i));
        end
        rhs(i) = F * a(i,:) * zj';
        
        % the resitive currents leaving the compartment
        jr(i) = sum((v(i)-v) ./ r(:,i));
    end
    
    % total current leaving this node from active and resistive processes.
    rhs = rhs + jr;

    dvdt = [0; dvdt_inv * rhs(2:end)];
end
