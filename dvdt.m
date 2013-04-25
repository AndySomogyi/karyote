function [dvdt] = dvdt(a, z, j, cap)
% An eqn for each compartment, a is define by
% F \sum_{a'!=a} A_{a,a'} z \cdot J_{a,a'} = 
% \sum_{a'!=a} A_{a,a'}C_{a,a'} \frac{dV_{a,a'}{dt}
% The dv/dt terms can be solved by treating this as a linear system, 
% Ax=b, where x is the column vector of dv/dt and solving for it. 
%
% a: area of membrane separating compartment a from a'
% z: valence of species i
% j: intra compartment flux, n_species,n_comp,n_comp
%
% returns a row vector of dvdt
%
    global F
    n = size(a,1);
    
    % the matricies, a and c have zero diagonals
    % when worked out, we can see that each diagonal entry
    % in the combined matrix ac is simply the negative sum of all the 
    % other elements in that rowl
    ac = -a .* cap;
    for i = 1:n
        ac(i,i) = -ac(i,:) * ac(i,:)';
    end
    
    rhs = zeros(n,1);
    zj = zeros(1,n);
    for i = 1:n
        for k = 1:n
            zj(k) = z * j(:,i,k);
        end
        rhs(i) = F * a(i,:) * zj';
    end

    dvdt = -rhs \ ac;
end
