function [ jt ] = charge_neutral_flux(j, z )

    sz = size(j);
    disp(sz)

    jt = zeros(sz);

    for i1=1:sz(1)
       for i2=1:sz(2)
           jt(i1,i2,:) =  squeeze(j(i1,i2,:))' - squeeze(z(i1,i2,:))' .* ...
               (...
               (squeeze(z(i1,i2,:))' * squeeze(j(i1,i2,:))) / ...
               (squeeze(z(i1,i2,:))' * squeeze(z(i1,i2,:))) ...
               );
       end
    end
    
    % verify that z . j == 0
    
    for i1=1:sz(1)
       for i2=1:sz(2)
           disp(squeeze(z(i1,i2,:))' * squeeze(jt(i1,i2,:)))
       end
    end
end

