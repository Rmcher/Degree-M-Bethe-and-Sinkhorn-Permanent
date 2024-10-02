% Author: WU Binghong
% Date: 2024.Oct.02

% Normal Permanent in recursion
function perm = computeNormalPermanent(A)
    [n, ~] = size(A);
    if n == 1
        perm = A(1,1);
        return;
    end
    perm = 0;
    for j = 1:n
        minor = A(2:end, [1:j-1, j+1:end]);
        perm = perm + A(1,j) * computeNormalPermanent(minor);
    end
end