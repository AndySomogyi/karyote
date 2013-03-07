function [ jt ] = charge_neutral_flux(j, z )
% j: total flux, n_comp x n_comp x n_species
% z: is the vector of valences, dims: 1 x n_species

    sz = size(j);
    disp(sz)

    jt = zeros(sz);
    
    % squeeze((x(i,i,:)) gives a column vector

    for i1=1:sz(1)
       for i2=1:sz(2)
           zz = z * z';
           if zz ~= 0
               jj = squeeze(j(i1,i2,:))';
               jt(i1,i2,:) =  jj - z .* (z * jj')/(zz); 
           else
               jt(i1,i2,:) = 0;
           end
       end
    end
    
    %disp('jt: ')
    %disp(jt)
    
    % verify that z . j == 0
    
    %for i1=1:sz(1)
    %   for i2=1:sz(2)
    %       z * squeeze(jt(i1,i2,:))
    %   end
    %end
end

