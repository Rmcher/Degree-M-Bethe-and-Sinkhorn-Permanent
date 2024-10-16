function [permscSA, DAE, D, E, iteration] = computeSinkPermanent(A, n)
% computeSinkPermanent - Computes the permanent of matrix A after applying
% Sinkhorn decomposition.
    %
    % Syntax: [permscSA, DAE, D, E, iteration] = computeSinkPermanent(A, n)
    %
    % Inputs:
    %   A - The input matrix (n x n).
    %   n - The dimension of the matrix A (integer).
    %
    % Output:
    %   permSA - The computed permanent using perm(D) * perm(M) * perm(E).
    %   DAE - The double stochastic matrix after Sinkhorn decomposition.
    %   D - D matrix.
    %   E - E matrix.
    %   iteration - The time to converge.
    %
    % Author: WU Binghong
    % Date: 2024.Oct.15
    
    tolerance = 1e-6;
    maxDifference = Inf;
    iteration = 0;

    D = diag(ones(1, n));
    E = diag(ones(1, n));

    DAE = A;

    while maxDifference > tolerance
        iteration = iteration + 1;

        % Update by column
        colsum = sum(DAE, 1);
        DAE = DAE ./ repmat(colsum, size(DAE, 1), 1);

        E = diag(colsum) * E;

        % Update by row
        rowsum = sum(DAE, 2);
        DAE = DAE ./ repmat(rowsum, 1, size(DAE, 2));

        D = D * diag(rowsum);

        % Check convergence
        maxDifference = max(abs(colsum - 1)) + max(abs(rowsum - 1));
    end

    % Compute permanents of D, M, and E
    permD = prod(diag(D));
    permE = prod(diag(E));

    % Compute final result
    permscSA = exp(-n) * permD * permE;
    
end