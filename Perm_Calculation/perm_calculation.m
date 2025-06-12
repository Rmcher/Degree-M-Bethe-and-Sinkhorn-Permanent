% Author: WU Binghong
% Date: 2024.Oct.02

% Define the number of matrices
numMatrices = 10;
n = 4;
M = 2;

% Define the arrays to store results
permValues = zeros(1, numMatrices);
permBValues = zeros(1, numMatrices);
permBMValues = zeros(1, numMatrices);
permscSValues = zeros(1, numMatrices);
permscSMValues = zeros(1, numMatrices);

% Generate numMatrices random matrices and calculate relevant results
for i = 1:numMatrices
    % Generate a n*n matrix with entries from uniform distribution [0, 1]
    A = rand(n, n);
    
    % Compute the permanent of A
    permValues(i) = perm(A, n);
    
    % Compute Bethe permanent approximations of A
    permBValues(i) = computeBethePermanent(A, n);

    % Compute Sinkhorn permanent approximations of A
    permscSValues(i) = computeSinkPermanent(A, n);

    % Compute M cover Bethe permanent approximations of A
    permBMValues(i) = computeBethePermanentM(A, n, M);

    % Compute M cover Sinkhorn permanent approximations of A
    permscSMValues(i) = computeSinkPermanentM(A, n, M);

    disp(i);
end

% Plot the results
figure;
hold on;

% Plot perm(A) vs perm_{B,M=2}(A) (red circles)
loglog(permValues, permBMValues, 'ro', 'DisplayName', sprintf('perm_{B, M=%d}(A)', M));

% Plot perm(A) vs perm_{S,M=2}(A) (yellow diamond)
loglog(permValues, permscSMValues, 'yd', 'DisplayName', sprintf('perm_{scS, M=%d}(A)', M));

% Plot perm(A) vs perm_B(A) (cyan triangles)
loglog(permValues, permBValues, 'c^', 'DisplayName', 'perm_B(A)');

% Plot perm(A) vs perm_S(A) (blue plus)
loglog(permValues, permscSValues, 'b+', 'DisplayName', 'perm_{scS}(A)');

% Plot
x_log = log10(permValues);
y_log = log10(permBValues);

x_center = median(x_log);
y_center = median(y_log);
W = 1.5;

theory_ratio = sqrt(exp(1) / (pi * n));
x_vals = logspace(x_center - W, x_center + W, 500);
y_vals = theory_ratio * x_vals;
loglog(x_vals, y_vals, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('y = %.3f x (theory 1)', theory_ratio));

xlim([10^(x_center - W), 10^(x_center + W)]);
ylim([10^(y_center - W), 10^(y_center + W)]);

% Set labels and legend
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('perm(A)');
ylabel('perm_B(A), perm_{B,M}(A), perm_{scS}(A) and perm_{scS,M}(A)');
legend;
title('Comparison of perm(A) vs perm_B(A), perm_{B,M}(A), perm_{scS}(A) and perm_{scS,M}(A)');

hold off;

% Save the plot as a .fig file with a name reflecting n and M
filename = sprintf('results/perm_comparison_n_%d_M_%d_num_%d.fig', n, M, numMatrices);
savefig(filename);
