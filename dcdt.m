function [ dcdt ] = dcdt( o, a, jt, nu, w )
%UNTITLED Summary of this function goes here
%    jt, charge-neutral flux of species from a' to a, dims: 
%       ncomp x ncomp x nxpecies
%    a: area of membrane separating a from a', dims: ncomp x ncomp
%    o: volume of compartments, dims: ncomp
%    w: rate of reaction k in compartment z, dims: n_comp x n_reactions
%    nu: = stoichiometric coefficient of species i for reaction k in
%        compartment a. dims: nreactions x ncomp x nspecies

    nspecies = size(jt,3);
    ncomp = size(jt,1);
    dcdt = zeros(ncomp,nspecies);
    
    for i = 1:ncomp
        %aa = a(i,:)               % row vector  (1, ncomp)
        %jj = squeeze(jt(i,:,:))   % matrix         (n_comp, n_species)
        %aj = 1/o(i) * a(i,:) * squeeze(jt(i,:,:)); %(1, n_species)
        dcdt(i,:) = 1/o(i) * a(i,:) * squeeze(jt(i,:,:));
        
        wk = w(i,:);                % row vector, (1, nreactions)
        nuik = squeeze(nu(:,i,:));  % matrix, n_reactions x n_species
        
        %disp(nuik)
        %disp(wk)
        %disp(wk * nuik)
        dcdt(i,:) = dcdt(i,:) + wk * nuik;
    end
end

