function [ y ] = affinity(nu, mu)
% y_{k,a} = \sum_i \nu_{ika} \mu_{ia}
% nu: stoichiometric cooefecients, i species x k reactions x a compartments
% mu electrochemical potential, i species x a compartments

    nu_size = size(nu);
    y_size = [nu_size(2), nu_size(3)];
    y = zeros(y_size);

    for k=1:y_size(1)
        for a = 1:y_size(2)
            y(k,a) = squeeze(nu(:,k,a))' * mu(:,a);
        end
    end
end

