function [ c, v ] = karyote_unpack(packed, n_comp, n_species)
% Packs the concentration and voltage matricies into a single column
% vector. 

% packed storage: first n_comp*n_species is the concentration matrix, 
% last n_comp is the compartment voltage. 

    n_time = size(packed, 1);
    c = zeros([n_time, n_species, n_comp]);
    % the packed time series of concentrations
    cpacked = packed(:, 1:(n_comp*n_species));
    
    for t=1:n_time
        cc = cpacked(t,:);
        c(t,:,:) = reshape(cc, [n_species, n_comp]);
    end

   
    v = packed(:, (n_comp*n_species)+1:end);


end

