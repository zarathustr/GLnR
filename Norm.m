function n = Norm(a, b)
s = size(a);
n = (trace(a' * a))^(1 / 2);
end