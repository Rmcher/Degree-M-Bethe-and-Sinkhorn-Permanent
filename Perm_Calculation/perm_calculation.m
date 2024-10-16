% Author: WU Binghong
% Date: 2024.Oct.02

% Define the number of matrices
numMatrices = 2;
n = 3;
M = 2;
% Define the arrays to store results
permValues = zeros(1, numMatrices);
permBValues = zeros(1, numMatrices);
permBMValues = zeros(1, numMatrices);
permSValues = zeros(1, numMatrices);
permSMValues = zeros(1, numMatrices);

% Generate 1000 random matrices and calculate relevant results
for i = 1:numMatrices
    % Generate a n*n matrix with entries from uniform distribution [0, 1]
    A = rand(n, n);
    
    % Compute the permanent of A
    permValues(i) = perm(A, n);
    
    % Compute Bethe permanent approximations of A
    permBValues(i) = computeBethePermanent(A, n);

    % Compute Sinkhorn permanent approximations of A
    permSValues(i) = computeSinkPermanent(A, n);

    % Compute M cover Bethe permanent approximations of A
    permBMValues(i) = computeBethePermanentM(A, n, M);

    % Compute M cover Sinkhorn permanent approximations of A
    permSMValues(i) = computeSinkPermanentM(A, n, M);
end

% Plot the results
figure;
hold on;

% Plot perm(A) vs perm_B(A) (cyan triangles)
loglog(permValues, permBValues, 'c^', 'DisplayName', 'perm_B(A)');

% Plot perm(A) vs perm_S(A) (blue plus)
loglog(permValues, permSValues, 'b+', 'DisplayName', 'perm_S(A)');

% Plot perm(A) vs perm_{B,M=2}(A) (red circles)
loglog(permValues, permBMValues, 'ro', 'DisplayName', sprintf('perm_{B, M=%d}(A)', M));

% Plot perm(A) vs perm_{S,M=2}(A) (yellow diamond)
loglog(permValues, permSMValues, 'yd', 'DisplayName', sprintf('perm_{S, M=%d}(A)', M));

% Set labels and legend
xlim([10^-(1.3) 10^(1.3)]);
ylim([10^-(1.3) 10^(1.3)]);
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('perm(A)');
ylabel('perm_B(A), perm_{B,M}(A), perm_S(A) and perm_{S,M}(A)');
legend;
title('Comparison of perm(A) vs perm_B(A), perm_{B,M}(A), perm_S(A) and perm_{S,M}(A)');

hold off;

% Save the plot as a .fig file with a name reflecting n and M
filename = sprintf('results/perm_comparison_n_%d_M_%d_num_%d.fig', n, M, numMatrices);
savefig(filename);