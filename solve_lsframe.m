function [W_store, M_store] = solve_lsframe(E_cov, K_cov)
    F = E_cov;
    [V_K, D_K] = eig(K_cov);
    [~,I_K] = sort(diag(D_K),'descend');
    V_K = V_K(:,I_K);
    D_K_diag = diag(D_K);
    D_K = diag(D_K_diag(I_K));
    R = V_K*inv(sqrt(D_K))*V_K.';
    M_store = {};
    W_store = {};
    % M init
    [M_old,D_M_old] = eig(R.'*F*R);
    [~,I_M_old] = sort(diag(D_M_old),'descend');
    M_store{1} = M_old(:,I_M_old);
    W_store{1} = nan;
    % Iteration
    M_diff = inf;
    iter_n = 1;
    while 1
        iter_n = iter_n + 1;
        if M_diff<1e-10 || iter_n>=100+1
            break;
        end
        % Update W
        W = inv(F) * F * R * M_store{iter_n-1};
        W_store{iter_n} = W;
        % Update M
        P = R.' * F * W;
        [U_P, ~, V_P] = svd(P,"econ");
        M = U_P * V_P.';
        M_store{iter_n} = M;
        % Update M_diff
        M_diff = mean(mean((M_store{iter_n} - M_store{iter_n-1}).^2));
    end
end