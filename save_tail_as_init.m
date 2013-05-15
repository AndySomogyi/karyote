function [  ] = save_tail_as_init( fname, c, v )
% Saves the last frame out of c and v as the variables
% 'c0' and 'v0' in a file with the given filename. 


    n_comp = size(v,2);
    n_species = size(c,2) / n_comp;
    
    clast = c(end,:);
    v0 = v(end,:);
    
    c0 = reshape(clast, [n_species, n_comp]);
    
    save(fname, 'c0', 'v0');
    
end

