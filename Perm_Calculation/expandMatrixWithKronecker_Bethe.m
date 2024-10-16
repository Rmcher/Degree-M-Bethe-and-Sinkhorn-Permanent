function expandedMatrices = expandMatrixWithKronecker_Bethe(A, n, M)
% expandMatrixWithKronecker_Bethe - Expands matrix A to dimension nM * nM
% for M cover Bethe approximation
    %
    % Syntax: expandedMatrices = expandMatrixWithKronecker_Bethe(A, n, M)
    %
    % Inputs:
    %   A - The input matrix (n x n).
    %   n - The dimension of the matrix A (integer).
    %   M - The expand dimension, cover of graphs.
    % Output:
    %   expandedMatrices - all possible M cover in cell.
    %
    % Author: WU Binghong
    % Date: 2024.Oct.15

    permutations = perms(1:M);
    numSpecialMatrices = size(permutations, 1);
    specialMatrices = cell(numSpecialMatrices, 1);

    for k = 1:numSpecialMatrices
        P = zeros(M);
        for i = 1:M
            P(i, permutations(k, i)) = 1;
        end
        specialMatrices{k} = P;
    end

    % Calculate the total number of expanded matrices
    totalExpandedMatrices = numSpecialMatrices^((n-1)^2);
    expandedMatrices = cell(totalExpandedMatrices, 1);
    
    % Generate all combinations of special matrices for non-first row and non-first column
    combinationIndices = repmat({1:numSpecialMatrices}, 1, (n-1)^2);
    allCombinations = combvec(combinationIndices{:});
    numCombinations = size(allCombinations, 2);
    
    % Iterate through each combination to create expanded matrices
    for combinationIdx = 1:numCombinations
        combination = allCombinations(:, combinationIdx);
        expandedA = zeros(n * M);

        for i = 1:n
            for j = 1:n
                if i == 1 || j == 1
                    % For the first row and first column, do Kronecker product with identity
                    kronProduct = A(i, j) * eye(M);
                else
                    % Use the specific permutation matrix from the combination
                    kronIndex = (i-2) * (n-1) + (j-2) + 1;
                    kronProduct = A(i, j) * specialMatrices{combination(kronIndex)};
                end

                % Place Kronecker product outcome to the new matrix
                rowStart = (i - 1) * M + 1;
                rowEnd = i * M;
                colStart = (j - 1) * M + 1;
                colEnd = j * M;

                expandedA(rowStart:rowEnd, colStart:colEnd) = kronProduct;
            end
        end

        expandedMatrices{combinationIdx} = expandedA;
    end

    return;
end