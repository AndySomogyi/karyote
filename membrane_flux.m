function [ j ] = membrane_flux( h,e )
% Calculates the total membrane flux, j, size: n_species x n_comp x n_comp.
% h: membrane permeability. size: n_species x n_comp x n_comp.
% e: equilibrium factor. size: n_species x n_comp x n_comp.
% Flux = Flow rate / membrane area
% e: mol / (L * nm), h: nm/s j-> mol/(L s)

    j = h .* e;
end

