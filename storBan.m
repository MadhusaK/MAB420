function [banMatrix,cholBan,sol] = storBan(coMatrix,b)
%% Cholesky with Band Storage and RCM node rodering
% 1. Generates RCM node reordered coMatrix
% 2. Generates Banded Storage of reordered coMatrix
% 3. Computes Cholesky Factorization of Banded system
% 4. Computes solution from banded cholesky

% Paramters
sol = b';
n = 33;

%% RCM Node Re-ordering
permRCM = symrcm(coMatrix);

rcmReOrder = triu(coMatrix(permRCM,permRCM)); % Reordered coefficient matrix
% rcmReOrder = triu(coMatrix); % Reordered coefficient matrix
sol = b(permRCM); % Reordered b vector

figure(1)
spy(rcmReOrder)

%% Banded Storage

banMatrix = zeros(7,33);

for i = 1:7
    banMatrix(8-i,i:33) = diag(rcmReOrder,i-1)'
end

figure(2)
spy(banMatrix)

%% Cholesky Factorization

cholBan = banMatrix;

for i = 1:33
    
    if i < 27
        % j range for traversing the rows of coMatrix in banded form
        jRange = 7:-1:1;
    else
        % Upper bound of the banded matrix
        jRange = 7:-1:(i-26);
    end
    
    for j = jRange     
        if i ~= 1 && j ~= 1
            cholBan(j,i+max(jRange)-j) = cholBan(j,i+max(jRange)-j) - dot(cholBan(1:(max(jRange)-1),i),[zeros(max(jRange)-j,1);cholBan(1:j-1,i+max(jRange)-j)])
        end
        
        if j == 7
            cholDiv = sqrt(cholBan(7,i));
        end
        
        cholBan(j,i+max(jRange)-j) = cholBan(j,i+max(jRange)-j)/cholDiv;
    end
end

%% Forward Subsitution
for i = 1:n
    
    if i ~= 1 && i <= 7
        i
        sol(i) = sol(i) - dot(sol(1:(i-1)),cholBan(1:(i-1),i));
    end
    
    if i ~= 1 && i > 7
        sol((i-6):(i-1))
        cholBan(7:-1:2)
        sol(i) = sol(i) - dot(sol((i-6):(i-1)),cholBan(1:6,i));

    end
    sol(i) = sol(i)/cholBan(7,i);
end


%% Backward Subtitution

for i = n:-1:1
    if i > 27 && i ~=33
        itt = 1; %needed not to override i
        for j = 6:-1:(7-n+i)
            sol(i) = sol(i) - sol(i+itt)*cholBan(j,i+itt);
            itt = itt+1;
        end
    end
    
    if i <= 27
        itt = 1;
        for j = 6:-1:1
            sol(i) = sol(i) - sol(i+itt)*cholBan(j,i+itt);
            itt = itt+1;
        end
    end
        
    sol(i) = sol(i)/cholBan(7,i)
end

sol = sol(permRCM');

