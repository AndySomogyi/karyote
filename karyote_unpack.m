function [ c, v ] = karyote_unpack(packed, n_comp, n_species)
% Packs the concentration and voltage matricies into a single column
% vector. 

% packed storage: first n_comp*n_species is the concentration matrix, 
% last n_comp is the compartment voltage. 

    n_time = size(packed, 1);

    c = reshape(packed(:, 1:(n_comp*n_species)), [n_time, n_comp, n_species]);
    v = reshape(packed(:, (n_comp*n_species)+1:end), [n_time, 1, n_comp]);


end

