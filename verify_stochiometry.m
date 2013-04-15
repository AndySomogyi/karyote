function [  ] = verify_stochiometry( s, z )
% Verify the validity of the stoichiometric matrix, 
%   checks for mass and charge balance. 

    n_species = size(s, 1);
    n_comp = size(s,2);
    n_reactions = size(s,3);
    
    s = reshape(s, [n_species*n_comp, n_reactions]);
    
    for i=1:n_reactions
        assert(sum(s(:,i)) == 0, sprintf('mass balance not conserved in reaction %i', i));
    end
    disp('mass balance appears OK');
    
    
    
    assert(n_species == length(z), 'valences must be the same length as the number of species');
    
    z = kron(z, ones(1,n_comp));
    for i=1:n_reactions
        assert(z * s(:,i) == 0, sprintf('charge balance not conserved in reaction %i', i));
    end
    disp('charge balance appears OK');
    
    disp('stoichiometry matrix appears OK');
    
end

