function P = P_matrix(x, ss)
    xx = x;
    P = zeros(ss, ss * (ss - 1) / 2);
    last_p = 1;
    for i = 1 : ss
        s = length(xx);
        p = last_p + s - 1;
        M = zeros(ss, s - 1);
        M(ss - s + 1, :) = xx(1, 2 : s);
        M(ss - s + 2 : ss, 1 : s - 1) = - xx(1) * eye(s - 1, s - 1);
        P(1 : ss, last_p : p - 1) = M;
        xx = xx(1, 2 : s);
        last_p = p;
    end
end