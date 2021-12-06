clear all
close all
clc

warning ('on','all');

tic;
dim = 4;
len = dim * (dim - 1) / 2;

disp('MATLAB Code:');

str = 'syms';

g_str = 'g = [';
for i = 1 : len
    str1 = sprintf(" g%d ", i);
    str = strcat(str, str1);
    g_str = strcat(g_str, str1);
end
g_str = strcat(g_str, '];');

str = strcat(str, 'real');
disp(str);
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
disp(str);
step2 = str;

str = 'G = [';
disp(str);
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
    disp(ss);
end
str = strcat(str, '];');
disp('];');
step3 = str;

disp(newline);
disp('G + transpose(G) => ');

eval(step1)
eval(step2)
eval(step3)
G + transpose(G)

disp(newline);
disp('P => ');

eval(g_str);
eval(x_str);
eqns = [G * x == zeros(dim, 1)];
[P, b] = equationsToMatrix(eqns, g);
P


disp(newline);
disp('Verification: (P * g - G * x)^T => ');
transpose(P * g' - G * x)
disp(newline);
generation_time = toc;
disp(sprintf('Generation time: %f s', generation_time));
disp(newline);

BB = randn(dim, dim);
[U, S, V] = svd(BB);
D = eye(dim, dim);
D(dim, dim) = det(U) * det(V);
C = V' * D * U;

B = zeros(dim, dim);
time_SVD = 0;
time_pro = 0;

num = floor(dim * 2);
H = zeros(len, len);
v = zeros(len, 1);
Db = zeros(num, dim);
Dr = randn(num, dim);
% Dr = ones(num, dim);
weights = abs(randn(num, 1));
weights = weights ./ sum(weights);
AA = [];
AAT = [];
BB = [];

U1 = [];
U2 = [];
U3 = [];
U4 = [];

scale = 1e-5;
Sigma_p = scale^2 * eye(dim);

T = randn(dim, 1) * 5;

mean_b = zeros(dim, 1);
mean_r = zeros(dim, 1);

for i = 1 : num
    r = Dr(i, :)';
    
    b = C * r + scale * randn(dim, 1) + T;
    
    Db(i, :) = b';
    Dr(i, :) = r';

    mean_b = mean_b + weights(i) * b;
    mean_r = mean_r + weights(i) * r;
end

for i = 1 : num
    b = Db(i, :)' - mean_b;
    r = Dr(i, :)' - mean_r;

    tic;
    s = b + r;
    d = b - r;
    PP = P_matrix(s', dim);
%     for j = 1 : dim
%         eval(sprintf('x%d = s(%d);', j, j));
%     end
%     PP = double(vpa(subs(P), 40));

    H = H + weights(i) * PP' * PP;
    v = v + weights(i) * PP' * d;
    AA = [AA, sqrt(weights(i)) * PP'];
    AAT = [AAT; sqrt(weights(i)) * PP];
    BB = [BB; sqrt(weights(i)) * d];
%     U1 = [U1; sqrt(weights(i)) * s(1) * eye(dim)];
%     U2 = [U2; sqrt(weights(i)) * s(2) * eye(dim)];
%     U3 = [U3; sqrt(weights(i)) * s(3) * eye(dim)];
%     U4 = [U4; sqrt(weights(i)) * s(4) * eye(dim)];
    time_pro = time_pro + toc;
    
    tic;
    B = B + weights(i) * b * r';
    time_SVD = time_SVD + toc;
end

% disp(sprintf('Determinant of H: %e', det(H)));

opts.tol = 1e-30;
opts.maxit = 5000;

tic;
for i = 1 : 10
    [U, S, V] = svds(B, dim, 'L', opts);
%     [U, S, V] = jacobi_svd(B);
    D = eye(dim, dim);
    D(dim, dim) = det(U) * det(V);
    CC = U * D * V';
end
disp(sprintf('MATLAB SVD Calculation Time: %5.8f s', toc));

% tic;
% for i = 1 : 1000
%     R_LMI = LMI(Db, Dr, weights);
% end
% disp(sprintf('LMI Calculation Time: %5.8f s', toc));

tic;
for i = 1 : 10
    tic;
    gg = H \ v;
end
disp(sprintf('Proposed Calculation Time: %5.8f s', toc));

for j = 1 : len
    eval(sprintf('g%d = gg(%d);', j, j));
end
GG = double(vpa(subs(G), 15));
C_res = (eye(dim) + GG) * ((eye(dim) - GG) \ eye(dim));

diff = norm(C - CC, 'inf');
disp(sprintf('SVD Difference between reference and restored rotation matrices: %5.16e', diff));

% diff = norm(C - R_LMI, 'inf');
% disp(sprintf('LMI Difference between reference and restored rotation matrices: %5.16e', diff));

diff = norm(C - C_res, 'inf');
disp(sprintf('Proposed Difference between reference and restored rotation matrices: %5.16e', diff));


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
disp(str);
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
variance_g = sqrt(diag(Sigma_g))';
R = C_res + eye(dim, dim);
K = zeros(dim, dim);
for i = 1 : dim
    sk = R(:, i);
    for j = 1 : dim
        str1 = sprintf(" x%d  = sk(%d);", j, j);
        eval(str1);
    end
    PD = eval(P);
    K = K + PD * Sigma_g * PD';
end

GGG = inv(eye(dim) + GG);
Sigma_R = GGG * K * GGG';

log(trace(Sigma_R)) / log(10.0)

% 
% A = H(1 : 3, 1 : 3);
% B = H(1 : 3, 4 : 6);
% C = H(4 : 6, 1 : 3);
% D = H(4 : 6, 4 : 6);
% 
% inv(A - B * inv(D) * C)
% inv(D) * C * inv(A - B * inv(D) * C)
% 
% U = chol(H);
% invU = inv(U);
% norm(inv(H) - invU * invU')
% 
% norm(pinv(BB) - BB' / norm(BB)^2)
% 
% norm(inv(AA * AAT)' * v - gg)
% 
% norm(pinv(AAT) - AA * pinv(AAT * AA))
% 
% norm(pinv(AA) - pinv(AAT * AA) * AAT)
% 
% norm(pinv(AA) - inv(AAT * AA)' * AAT)
% 
% norm(pinv(AAT) - pinv(AA)')

% simplify((-P'*P+(x'*x-1)*eye(len))/(x'*x-1) - simplify(pinv(eye(len) - P' * P)))

% P1 = [0, 0, 0, 0, 0, 0;
%       -1, 0, 0, 0, 0, 0;
%       0, -1, 0, 0, 0, 0;
%       0, 0, -1, 0, 0, 0];
%   
% P2 = [1, 0, 0, 0, 0, 0;
%       0, 0, 0, 0, 0, 0;
%       0, 0, 0, -1, 0, 0;
%       0, 0, 0, 0, -1, 0];
%   
% P3 = [0, 1, 0, 0, 0, 0;
%       0, 0, 0, 1, 0, 0;
%       0, 0, 0, 0, 0, 0;
%       0, 0, 0, 0, 0, -1];
%   
% P4 = [0, 0, 1, 0, 0, 0;
%       0, 0, 0, 0, 1, 0;
%       0, 0, 0, 0, 0, 1;
%       0, 0, 0, 0, 0, 0];
%   
%   
% BB' / norm(BB)^2 * (U1 * P1 + U2 * P2 + U3 * P3 + U4 * P4)
% gg'


% W = P' * P; W = W(1 : dim - 1, 1 : dim - 1); simplify(inv(W) * det(W))
% W
% 
% W = P' * P; W = W(dim : len, dim : len); simplify(inv(W) * det(W))
% W