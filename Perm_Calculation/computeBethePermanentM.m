function permBM = computeBethePermanentM(A, n, M)
% computeBethePermanentM - compute the M cover Bethe approximation
% permanent of matrix A based on standard permanent calculation
    %
    % Syntax: permBM = computeBethePermanentM(A, n, M)
    %
    % Inputs:
    %   A - The input matrix (n x n).
    %   n - The dimension of the matrix A (integer).
    %   M - Expand dimension, M cover of graphs.
    %
    % Output:
    %   permBM - The computed permanent using M cover Bethe Approximation.
    %
    % Author: WU Binghong
    % Date: 2024.Oct.15


    expandedMatrices = expandMatrixWithKronecker_Bethe(A, n, M);
    numExpandedMatrices = length(expandedMatrices);

    % Initialize vector to store permB of expanded matrices
    permBValues = zeros(numExpandedMatrices, 1);

    n_expanded = n * M;

    % Calculate permB for each matrix
    for idx = 1:numExpandedMatrices
        % Get each expanded matrix
        A_expanded = expandedMatrices{idx};

        % Compute permB in standard way
        permB = perm(A_expanded, n_expanded);

        % Store permB
        permBValues(idx) = permB;
    
        % Display progress
        % fprintf('%d of %d perm completed.\n', idx, numExpandedMatrices);
    end


    % Calculate permBM
    % sumperm = sum(permBValues);
    permBM = ((1 / (factorial(M)^((n-1)^2))) * sum(permBValues)) ^ (1 / M);

    return;
end
