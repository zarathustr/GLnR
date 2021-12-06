% Generalized Linear n-Dimensional Registration Algorithm
% By Jin Wu
% (c) 2018 www.jinwu.science
% jin_wu_uestc@hotmail.com


function [G, P, g, x] = GP_generator(dim)
len = dim * (dim - 1) / 2;

str = 'syms';

g_str = 'g = [';
for i = 1 : len
    str1 = sprintf(" g%d ", i);
    str = strcat(str, str1);
    g_str = strcat(g_str, str1);
end
g_str = strcat(g_str, '];');

str = strcat(str, 'real');
step1 = str;

str = 'syms';

x_str = 'x = [';
for i = 1 : dim
    str1 = sprintf(" x%d ", i);
    str = strcat(str, str1);
    x_str = strcat(x_str, str1);
end
x_str = strcat(x_str, ']'';');

str = strcat(str, 'real');
step2 = str;

str = 'G = [';
for i = 1 : dim
    ss = '';
    for j = 1 : dim
        if(i == j)
            str1 = sprintf(" 0");
        elseif(i < j)
            str1 = sprintf(" + g%d", (i - 1) * dim - i * (i + 1) / 2 + j);
        else
            str1 = sprintf(" - g%d", (j - 1) * dim - j * (j + 1) / 2 + i);
        end
        
        str = strcat(str, str1);
        ss = strcat(ss, str1);
        if(j < dim)
            str = strcat(str, ',');
            ss = strcat(ss, ',');
        end
    end
    str = strcat(str, ';');
    ss = strcat(ss, ';');
end
str = strcat(str, '];');
step3 = str;

eval(step1);
eval(step2);
eval(step3);

eval(g_str);
eval(x_str);
eqns = [G * x == zeros(dim, 1)];
[P, b] = equationsToMatrix(eqns, g);
end