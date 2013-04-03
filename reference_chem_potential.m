function [ mu ] = reference_chem_potential( k )
% Calculate the reference chemical potential based on the
% rate constants.  
% 
% k: forward and backward rate constants, n_reactions x 2
% returns: n_reaction x 1 vector, one for each reaction. 

    global R T;
    
    mu =  R * T * log(k(:,1) ./ k(:,2));
end

