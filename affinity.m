function [ y ] = affinity(nu, mu)
% y_{k,a} = \sum_i \nu_{ika} \mu_{ia}
% nu: stoichiometric cooefecients, nu = (n_reactions, n_comp, n_species)
% mu electrochemical potential, n_comp x n_species
% return, y: n_comp x n_reactions

    n_comp = size(mu, 1);
    y = zeros(n_comp, size(nu,1));
    
    for i = 1:n_comp
        vik = squeeze(nu(:,i,:))'; % n_species x n_reactions
        mui = mu(i,:);             % 1 x n_species
        y(i,:) = mui * vik;
    end
end

