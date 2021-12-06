function [R, T] = LMI(Db, Dr, weights)
    s = size(Db);
    num = s(1);
    dim = s(2);
    mean_b = zeros(dim, 1);
    mean_r = zeros(dim, 1);

    for i = 1 : num
        r = Dr(i, :)';
        b = Db(i, :)';

        mean_b = mean_b + weights(i) * b;
        mean_r = mean_r + weights(i) * r;
    end
    
    B = zeros(dim, dim);
    for i = 1 : num
        b = Db(i, :)' - mean_b;
        r = Dr(i, :)' - mean_r;
        B = B + weights(i) * b * r';
    end
    
    setlmis([]); 
    C = lmivar(2, [dim dim]);
    
    lmiterm([-1 1 1 0], eye(dim))
    lmiterm([-1 1 2 -C], eye(dim), eye(dim))
    lmiterm([-1 2 1 C], eye(dim), eye(dim))
    lmiterm([-1 2 2 0], eye(dim))
    LMIs = getlmis;
    c = mat2dec(LMIs, - B);
    opt = [1e-5, 0, -1, 0, 1];
    [copt, xopt] = mincx(LMIs, c, opt);
    R = reshape(xopt, dim, dim)';
    [u, ~, v] = svd(R);
    d = eye(dim);
    d(dim, dim) = det(u * v);
    R = u * d * v';
end