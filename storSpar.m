function [sparMatrix,cholSpar,solSpar] = storSpar(coMatrix,b)
%% Cholesky with Sparse Storage and AMD node rodering
% 1. Generates AMD node reordered coMatrix
% 2. Generates Sparse Storage of reordered coMatrix
% 3. Computes Cholesky Factorization of Sparse system
% 4. Computes solution from banded cholesky
%% Parameters
n = 33;
solSpar = b;
coMatrix = triu(coMatrix);
%% AMD Node reordering

permAMD = symamd(coMatrix);

amdReOrder = triu(coMatrix(permAMD,permAMD));

spy(amdReOrder)

sparMatrix = sparse(amdReOrder);

cholSpar = sparMatrix;
tic;
for i = 1:n
    if i ~= 1
        cholSpar(i,i:n) = cholSpar(i,i:n) - cholSpar(1:i-1,i)'*cholSpar(1:i-1,i:n);
    end
    cholSpar(i,i:n) = cholSpar(i,i:n)/sqrt(cholSpar(i,i));
end
toc;

% cholSpar = full(cholSpar);


%% Forward Substituion Gx = b
for i = 1:n
    if i ~= 1
        solSpar(i) = solSpar(i) - dot(cholSpar(1:i-1,i),solSpar(1:i-1));
    end
    
    solSpar(i) = solSpar(i)/cholSpar(i,i);
end


% 
% %% Back Subsitution G'z=x
% 
for i = n:-1:1
    if i ~= n
        solSpar(i) = solSpar(i) - dot(cholSpar(i,i+1:n),solSpar(i+1:n));
    end
    
    solSpar(i) = solSpar(i)/cholSpar(i,i);
end
