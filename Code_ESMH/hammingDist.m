
function Dh=hammingDist(B1, B2)  %%数据库、查询图像的哈希码   一行为一个哈希码
%%panpan
% A矩阵的每一行分别与B矩阵每行计算汉明距离  然后结果排列成矩阵
% Compute hamming distance between two sets of samples (B1, B2)
%
% Dh=hammingDist(B1, B2);
%
% Input
%    B1, B2: compact bit vectors. Each datapoint is one row.
%    size(B1) = [ndatapoints1, nwords]
%    size(B2) = [ndatapoints2, nwords]
%    It is faster if ndatapoints1 < ndatapoints2
% 
% Output
%    Dh = hamming distance. 
%    size(Dh) = [ndatapoints1, ndatapoints2]

% example query
% Dhamm = hammingDist(B2, B1);
% this will give the same result than:
%    Dhamm = distMat(U2>0, U1>0).^2;
% the size of the distance matrix is:
%    size(Dhamm) = [Ntest x Ntraining]

% loop-up table:

bit_in_char = uint16([...
    0 1 1 2 1 2 2 3 1 2 2 3 2 3 3 4 1 2 2 3 2 3 ...
    3 4 2 3 3 4 3 4 4 5 1 2 2 3 2 3 3 4 2 3 3 4 ...
    3 4 4 5 2 3 3 4 3 4 4 5 3 4 4 5 4 5 5 6 1 2 ...
    2 3 2 3 3 4 2 3 3 4 3 4 4 5 2 3 3 4 3 4 4 5 ...
    3 4 4 5 4 5 5 6 2 3 3 4 3 4 4 5 3 4 4 5 4 5 ...
    5 6 3 4 4 5 4 5 5 6 4 5 5 6 5 6 6 7 1 2 2 3 ...
    2 3 3 4 2 3 3 4 3 4 4 5 2 3 3 4 3 4 4 5 3 4 ...
    4 5 4 5 5 6 2 3 3 4 3 4 4 5 3 4 4 5 4 5 5 6 ...
    3 4 4 5 4 5 5 6 4 5 5 6 5 6 6 7 2 3 3 4 3 4 ...
    4 5 3 4 4 5 4 5 5 6 3 4 4 5 4 5 5 6 4 5 5 6 ...
    5 6 6 7 3 4 4 5 4 5 5 6 4 5 5 6 5 6 6 7 4 5 ...
    5 6 5 6 6 7 5 6 6 7 6 7 7 8]);

n1 = size(B1,1);
[n2, nwords] = size(B2); %n1为数据库的样本个数    n2为查询的样本个数

Dh = zeros([n1 n2], 'uint16');
for j = 1:n1            %n1个样本个数
    for n=1:nwords      %哈希码位数
        y = bitxor(B1(j,n),B2(:,n));    %
        Dh(j,:) = Dh(j,:) + bit_in_char(y+1);
    end
end
end
