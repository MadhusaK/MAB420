function [sparMatrix,cholSpar] = storSpar(coMatrix)
%% Cholesky with Sparse Storage and AMD node rodering
% 1. Generates AMD node reordered coMatrix
% 2. Generates Sparse Storage of reordered coMatrix
% 3. Computes Cholesky Factorization of Sparse system
% 4. Computes solution from banded cholesky
%% Parameters
n = 33;
%% AMD Node reordering

permAMD = symamd(coMatrix);

amdReOrder = triu(coMatrix(permAMD,permAMD));

spy(amdReOrder)

sparMatrix = sparse(amdReOrder);

cholSpar = sparMatrix;
for i = 1:n
    if i ~= 1
        cholSpar(i,i:n) = (cholSpar(i,i:n) - cholSpar(1:i-1,i)'*cholSpar(1:i-1,i:n));
    end
    cholSpar(i,i:n) = cholSpar(i,i:n)/sqrt(cholSpar(i,i));
end

