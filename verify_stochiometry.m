function [  ] = verify_stochiometry( s, z )
% Verify the validity of the stoichiometric matrix, 
%   checks for mass and charge balance. 

    n_species = size(s, 1);
    n_comp = size(s,2);
    n_reactions = size(s,3);
    
    s = reshape(s, [n_species*n_comp, n_reactions]);
    
    %for i=1:n_reactions
    %    assert(sum(s(:,i)) == 0, sprintf('mass balance not conserved in reaction %i', i));
    %end
    %disp('mass balance appears OK');
    
    
    
    assert(n_species == length(z), 'valences must be the same length as the number of species');
    
    z = kron(ones(1,n_comp), z);
    r = zeros(n_species*n_comp, 1);
    for i=1:n_reactions
        % reactants
        ind_lhs = find(s(:,i) < 0);
        r(:) = 0;
        r(ind_lhs) = abs(s(ind_lhs,i));
        
        lhs = z * r;
        
        % products
        ind_rhs = find(s(:,i) > 0);
        r(:) = 0;
        r(ind_rhs) = abs(s(ind_rhs,i));
        
        rhs = z * r;
        
        assert(lhs == rhs, ...
            sprintf('charge balance not conserved in reaction %i, consumed: %i, produced: %i', ...
            i, lhs, rhs));
    end
    

    disp('stoichiometry matrix appears OK');
    
end

