clear all
close all
clc

warning ('off','all');
Start = 3;
End = 200;
tt = zeros(End - Start + 1, 1);

for kk = Start : End
    tic;
    dim = 3;
    len = dim * (dim - 1) / 2;

    % disp('MATLAB Code:');

    str = 'syms';

    g_str = 'g = [';
    for i = 1 : len
        str1 = sprintf(" g%d ", i);
        str = strcat(str, str1);
        g_str = strcat(g_str, str1);
    end
    g_str = strcat(g_str, '];');

    str = strcat(str, 'real');
    % disp(str)
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
    % disp(str)
    step2 = str;

    str = 'G = [';
    % disp(str);
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
        % disp(ss);
    end
    str = strcat(str, '];');
    % disp('];');
    step3 = str;

    % disp(newline);
    % disp('G + transpose(G) => ');

    eval(step1)
    eval(step2)
    eval(step3)
    % G + transpose(G)

    % disp(newline);
    % disp('P => ');

    eval(g_str);
    eval(x_str);
    eqns = [G * x == zeros(dim, 1)];
    [P, b] = equationsToMatrix(eqns, g);
    % P


    % disp(newline);
    % disp('Verification: (P * g - G * x)^T => ');
    % transpose(P * g' - G * x)
    % disp(newline);
    tt(kk - Start + 1) = toc;
    % disp(sprintf('Generation time: %f s', generation_time));
    % disp(newline);
end

figure(1);
plot(Start + 3 : End, tt(4 : length(tt), 1), '*-', 'LineWidth', 1);