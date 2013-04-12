function [ j ] = membrane_flux( h,e )
% Calculates the total membrane flux, j, size: n_species x n_comp x n_comp.
% h: membrane permeability. size: n_species x n_comp x n_comp.
% e: equilibrium factor. size: n_species x n_comp x n_comp.

    j = h .* e;
end

