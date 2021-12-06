% Generalized Linear n-Dimensional Registration Algorithm
% By Jin Wu
% (c) 2018 www.jinwu.science
% jin_wu_uestc@hotmail.com



function [R, T, B, metric_error, Sigma_g, Sigma_R, Sigma_T] = GLnR(Db, Dr, weights, Sigma_p)
    s = size(Db);
    num = s(1);
    dim = s(2);
    len = dim * (dim - 1) / 2;
    mean_b = zeros(dim, 1);
    mean_r = zeros(dim, 1);

    [G, P, g, x] = GP_generator(dim);

    for i = 1 : num
        r = Dr(i, :)';
        b = Db(i, :)';

        mean_b = mean_b + weights(i) * b;
        mean_r = mean_r + weights(i) * r;
    end
    
    H = zeros(len, len);
    v = zeros(len, 1);
    B = zeros(dim, dim);
    for i = 1 : num
        b = Db(i, :)' - mean_b;
        r = Dr(i, :)' - mean_r;
    
        s = b + r;
        d = b - r;
        for j = 1 : dim
            eval(sprintf('x%d = s(%d);', j, j));
        end
        PP = double(vpa(subs(P), 40));

        H = H + weights(i) * PP' * PP;
        v = v + weights(i) * PP' * d;
        B = B + weights(i) * b * r';
    end
    gg = H \ v;
    for j = 1 : len
        eval(sprintf('g%d = gg(%d);', j, j));
    end
    GG = double(vpa(subs(G), 15));
    R = (eye(dim) + GG) * inv(eye(dim) - GG);
    T = mean_b - R * mean_r;
    
    metric_error = 0;
    for i = 1 : num
        b = Db(i, :)';
        r = Dr(i, :)';
    
        metric_error = metric_error + weights(i) * norm(b - R * r - T)^2;
    end
    
    PX = P;
    str = 'syms';

    d_str = 'd = [';
    for i = 1 : dim
        str1 = sprintf(" d%d ", i);
        str = strcat(str, str1);
        d_str = strcat(d_str, str1);
    end
    d_str = strcat(d_str, ']'';');

    str = strcat(str, 'real');
    eval(str);
    eval(d_str);

    for i = 1 : dim
        str1 = sprintf(" x%d  = d%d;", i, i);
        eval(str1);
    end
    PD = eval(P);
    HH = zeros(len, len);

    for i = 1 : num
        b = Db(i, :)';
        r = Dr(i, :)';
    
        s = b + r;

        KK = zeros(len, dim);
        for j = 1 : len
            J = PX' * PD + PD' * PX;
            K = equationsToMatrix(J(:, j), d);
            for k = 1 : dim
                str1 = sprintf(" x%d  = s(%d);", k, k);
                eval(str1);
            end
            KK = KK + gg(j) * eval(K);
        end
        KK = weights(i) * KK;

        for k = 1 : dim
            str1 = sprintf(" x%d  = b(%d);", k, k);
            eval(str1);
        end
        PPP = eval(P);

        HH = HH + 4 * weights(i)^2 * PPP' * Sigma_p * PPP + KK * Sigma_p * KK';
    end

    Sigma_g = inv(H) * HH * inv(H);
    
    RR = R + eye(dim, dim);
    K = zeros(dim, dim);
    for i = 1 : dim
        sk = RR(:, i);
        for j = 1 : dim
            str1 = sprintf(" x%d  = sk(%d);", j, j);
            eval(str1);
        end
        PD = eval(P);
        K = K + PD * Sigma_g * PD';
    end

    GGG = inv(eye(dim) + GG);
    Sigma_R = GGG * K * GGG';
    
    K = zeros(dim, dim);
    for i = 1 : dim
        sk = RR(:, i);
        for j = 1 : dim
            str1 = sprintf(" x%d  = mean_r(%d) * sk(%d);", j, j, j);
            eval(str1);
        end
        PD = eval(P);
        K = K + PD * Sigma_g * PD';
    end
    Sigma_T = Sigma_p + R * Sigma_p * R' + GGG * K * GGG';
end