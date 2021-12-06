function [U, S, V] = jacobi_svd(A)
    TOL = 1e-16;
    [m, n] = size(A);
    U = A;
    V = eye(n);
    converge = TOL + 1;
    cnt = 0;
    
    while converge > TOL
        converge = 0;
        
        if(cnt > 60)
            break;
        end
    
        for j = 2 : n
            for i = 1 : j - 1
                % compute [alpha gamma;gamma beta]=(i,j) submatrix of U'*U
                alpha = U(:, i)' * U(:, i); %might be more than 1 line
                beta = U(:, j)' * U(:, j); %might be more than 1 line
                gamma = U(:, i)' * U(:, j); %might be more than 1 line
                converge = max(converge, abs(gamma));
                % compute Jacobi rotation that diagonalizes
                % [alpha gamma;gamma beta]
                if gamma == 0
                    c = 1;
                    s = 0;
                else
                    zeta = (beta - alpha) / (2 * gamma);
                    if norm(zeta) == 0
                        t = 1 / (abs(zeta) + sqrt(1 + zeta^2));
                    else
                        t = sign(zeta) / (abs(zeta) + sqrt(1 + zeta^2));
                    end
                    c = 1 / sqrt(1 + t^2);
                    s = c * t;
                end
                % update columns i and j of U
                t = U(:, i);
                U(:, i) = c * t - s * U(:, j);
                U(:, j) = s * t + c * U(:, j);
                % update matrix V of right singular vectors
                t = V(:, i);
                V(:, i) = c * t - s * V(:, j);
                V(:, j) = s * t + c * V(:, j);
            end
        end
        cnt = cnt + 1;
    end
    % the singular values are the norms of the columns of U
    % the left singular vectors are the normalized columns of U
    for j = 1 : n
        singvals(j) = norm(U(:, j));
        U(:, j) = U(:, j) / singvals(j);
    end
    [~, in] = sort(singvals, 2, 'descend');
    S = diag(singvals(in));
    V = V(:, in);
    r = rank(S);
    U = U(:, in);
    U1 = U(:, 1 : r);
    U2 = null(U1');
    U = [U1, U2];
    if m > n
        S(n + 1 : m, :) = zeros(m - n, n);
    end
end