function [pacVector,cholPac,sol] = storPac(coMatrix,b)
%% Cholesky with Packed Storage
% 1. Generates packstorage format for coMatrix
% 2. Computes the cholesky factorization for packed coMatrix
% 3. Computes solution from packed cholesky 

%% Paramters
cholPac = zeros(size(coMatrix));
n = 33;
upTri = triu(coMatrix);
pacVector = [];
sol = b;

%% Generate packed storage - column wise
for i = 1:33
        pacVector = [pacVector;upTri(1:i,i)];
end

%% Functoin: Coordinate -> Index

f_nCo = @(i,j) i+j*(j-1)/2;
f_nRo = @(i,j) j+i*(i-1)/2;

%% Cholesky Factorization
cholPac = pacVector;
 
tic;
for i = 1:n
    cholDiv = sqrt(cholPac(f_nCo(i,i)));
    for j = i:n
        if i ~= 1
            cholPac(f_nCo(i,j)) = cholPac(f_nCo(i,j)) - dot(cholPac(f_nCo(1,i):f_nCo(i-1,i)),cholPac(f_nCo(1,j):f_nCo(i-1,j)));
        end
        
        if i == j
            cholDiv = sqrt(cholPac(f_nCo(i,i)));
        end
        
        cholPac(f_nCo(i,j)) = cholPac(f_nCo(i,j))/cholDiv;
    end
end
time = toc;
fprintf('Packed Storage - Cholesky:             %10.10fms \n', time)
%% Forward Substitution
tic;
for i = 1:n
    for j = 1:i-1
        sol(i) = sol(i) - cholPac(f_nCo(j,i))*sol(j);
    end

    sol(i) = sol(i)/cholPac(f_nCo(i,i));
end
        
%% Backwards Substitution
for i = n:-1:1
    for j = i+1:n
        sol(i) = sol(i) - cholPac(f_nCo(i,j))*sol(j);
    end
    
    if isempty(j)
        j = 33;
    end
    
    sol(i) = sol(i)/cholPac(f_nCo(i,i));
end
time = toc;
fprintf('Packed Storage - Solution Runtime:     %10.10fms \n', time)

figure
spy(pacVector')
title('Packed Storage Vector')
