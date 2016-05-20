F = 50;
K = 20;
P = 40;

%% Sigle subspace
%figure;
%imshow(rand(F, K, 3));
%title('Valid Coefficient');
%
%figure;
%imshow(0.5 + 0.5*rand(F, P-K));
%title('invalid coefficient')
%
%figure;
%imshow(rand(K, P, 3));
%title('valid bases')
%
%figure;
%imshow(0.5 + 0.5*rand(P-K, P));
%title('invalid bases')
%
% Sparsity
%figure;
%L = 3*K;
%s = 5;
%im = 0.9+ 0.1*rand(F, L, 3);
%for i = 1:F
%    im(i, randsample(L, s), :) = rand(1, s, 3);
%end
%imshow(im);
%
%% Bases
%figure;
%imshow(rand(L, P, 3));
