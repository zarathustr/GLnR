function G = G_matrix(g, dim)
    G = zeros(dim, dim);
    last_p = 1;
    for i = 1 : dim - 2
        p = last_p + dim - i;
        G(i, i + 1 : dim) = g(last_p : p - 1);
        last_p = p;
    end
    G(dim - 1, dim) = g(dim * (dim - 1) / 2);
    
    for i = 1 : dim
        for j = 1 : dim
            G(j, i) = -G(i, j);
        end
    end
end