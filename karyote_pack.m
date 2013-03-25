function [ packed ] = karyote_pack( c,v )
% Packs the concentration and voltage matricies into a single column
% vector. 

% packed storage: first n_comp*n_species is the concentration matrix, 
% last n_comp is the compartment voltage. 

    n_comp = size(c,1);
    n_species = size(c,2);
    packed = zeros((n_comp*n_species) + n_comp, 1);

    
    %     c = reshape(s(1:(n_comp*n_species)), [n_comp, n_species]);
    % v = reshape(s((n_comp*n_species)+1:end), [1, n_comp]);

    packed(1:(n_comp*n_species)) = reshape(c, [1,(n_comp*n_species)]);
    packed((n_comp*n_species)+1:end) = v;


end

