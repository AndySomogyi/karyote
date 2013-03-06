function [ m ] = zero_diagonal( m )
%   zero the diagonal of any square matrix
    m(logical(eye(size(m)))) = 0;
end

