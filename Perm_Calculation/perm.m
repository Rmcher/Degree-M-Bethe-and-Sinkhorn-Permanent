function standardperm = perm(A, n)
% perm - Computes standard permanent base on math definition.
    %
    % Syntax: standardperm = perm(A)
    %
    % Inputs:
    %   A - The input matrix (n x n).
    %   n - The dimension of the matrix A (integer).
    %
    % Output:
    %   standardperm - The standard permanent.
    %
    % Author: WU Binghong
    % Date: 2024.Oct.02
    
    if n == 1
        standardperm = A(1,1);
        return;
    end
    standardperm = 0;
    for j = 1:n
        minor = A(2:end, [1:j-1, j+1:end]);
        standardperm = standardperm + A(1,j) * perm(minor, n-1);
    end
end