function [sol,cholFac,cholFac2] = storFull(coMatrix,b)
% Computes the cholesky facorization and solves the matrix

n = 33;
sol = b;
cholFac = coMatrix; % Non-vectorized Solution
cholFac2 = coMatrix; % Vectorized Solution

%% Vectorized Solution

tic;
for i = 1:n
    if i ~= 1
        cholFac2(i,i:n) = (cholFac2(i,i:n) - cholFac2(1:i-1,i)'*cholFac2(1:i-1,i:n));
    end
    cholFac2(i,i:n) = cholFac2(i,i:n)/sqrt(cholFac2(i,i));
end
time = toc;
fprintf('Full Storage - Vectorized Cholesky:    %10.10fms \n', time)


cholFac2 = triu(cholFac2); 
%% Cholesky Facorization
% Computes Cholesky with standard cholesky algorithm (slow -- well, if the matrix is large enough)

tic;
for j = 1:n
    for i = 1:j-1
        for k = 1:i-1
            cholFac(i,j) = cholFac(i,j) - cholFac(k,i)*cholFac(k,j);
        end
        
        cholFac(i,j) = cholFac(i,j)/cholFac(i,i);
        cholFac(j,j) = cholFac(j,j) - cholFac(i,j)^2;
        
    end
    cholFac(j,j) = sqrt(cholFac(j,j));
end
time = toc;
fprintf('Full Storage - Itt Cholesky:           %10.10fms \n', time)



%% Establish lower and uppper matrix
% cholFac = triu(cholFac);
cholFac = cholFac2;

tic;
%% Forward Substituion Gx = b
for i = 1:n
    if i ~= 1
        sol(i) = sol(i) - dot(cholFac(1:i-1,i),sol(1:i-1));
    end
    
    sol(i) = sol(i)/cholFac(i,i);
end


% 
% %% Back Subsitution G'z=x
% 
for i = n:-1:1
    if i ~= n
        sol(i) = sol(i) - dot(cholFac(i,i+1:n),sol(i+1:n));
    end
    
    sol(i) = sol(i)/cholFac(i,i);
end
time = toc;
fprintf('Full Storage - Solution Runtime:       %10.10fms \n', time)
