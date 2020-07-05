%% Test Givens
clc; clear;

M = zeros(2,4);
M(:) = 1:8;

for i=1:size(M,2)
    for j =size(M,1):-1:i+1
        [G,M] = givens(M,i,j);
    end
end

M

non_zero = find(~all(M<eps*1000,2))

M = M(non_zero,:)