function [coMatrix,coB] = genAB
% Constructs Coefficient matrix for Group Assignment


%Coef Vector:   U(0,0),...,U(0,5),  -- 6 Points
%               U(1,0),...,U(1,5),  -- 6 Points
%               U(2,0),...,U(2,5),  -- 6 Points
%               U(3,0),...,U(3,4),  -- 6 Points
%               U(4,0),...,U(4,3),  -- 5 Points
%               U(5,0),...,U(5,2),  -- 4 Points
%                                   -- 33 Points -> 33x33 Co.Matrix

% Initialize Coefficients and matricies
tInf = 20;
k = 3;
h = 20;
delX = 0.01;
alpha = (h*delX)/(k);
coMatrix = zeros(33,33);
coB = zeros(33,1);

%% Functoinl: Coordinate -> Vector Index  %%
fn_nPoints = @(i,j) i*6 + (j+1) - round(i/9)*(1);


%% Approximation for points U(i,0) : i = 0:5 %%

for i = 0:5
    coMatrix(fn_nPoints(i,0),fn_nPoints(i,0)) = 3 - floor(1/(i+1)); %U(i,0), floor(..) deals with the special case at U(0,0)
    coMatrix(fn_nPoints(i,0),fn_nPoints(i,1)) = -1; %U(i,1)
    
    if i<5 %U(i+1,0)
        coMatrix(fn_nPoints(i,0),fn_nPoints(i+1,0)) = -1; 
    end
    
    if i>0 %U(i-1,0)
        coMatrix(fn_nPoints(i,0),fn_nPoints(i-1,0)) = -1;
    end
end
    
%% Approximation for points U(i,j) : i = 0:5, j=1:2

for j = 1:2
    for i = 0:5
        coMatrix(fn_nPoints(i,j),fn_nPoints(i,j)) = 4 - floor(1/(i+1)); %U(i,j), floor(..) deals with the special case at U(0,0)
        coMatrix(fn_nPoints(i,j),fn_nPoints(i,j-1))= -1; %U(i,j-1)
        coMatrix(fn_nPoints(i,j),fn_nPoints(i,j+1)) = -1; %U(i,j+1)

        if i<5 %U(i+1,j)
            coMatrix(fn_nPoints(i,j),fn_nPoints(i+1,j)) = -1; 
        end

        if i>0 %(i-1,j)
            coMatrix(fn_nPoints(i,j),fn_nPoints(i-1,j)) = -1;
        end
    end
end
    
%% Approximation for points U(i,j): j = 3:5, i = 0:(7-j)

for j = 3:5
    for i = 0:(7-j)
        coMatrix(fn_nPoints(i,j),fn_nPoints(i,j)) = 4 - floor(1/(i+1)); %U(i,j), floor(..) deals with the special case at U(0,0)
        coMatrix(fn_nPoints(i,j),fn_nPoints(i,j-1))= -1 %U(i,j-1)
        
        if j<5 %U(i,j+1)
            coMatrix(fn_nPoints(i,j),fn_nPoints(i,j+1)) = -1;         
        end
        
        if i<5 %U(i+1,j)
            coMatrix(fn_nPoints(i,j),fn_nPoints(i+1,j)) = -1; 
        end        
        
        if i>0 %(i-1,j)
            coMatrix(fn_nPoints(i,j),fn_nPoints(i-1,j)) = -1;
        end
    end
end

%% Approximation for U(5,3), U(4,4), and U(3,5)

%U(5,3)
coMatrix(fn_nPoints(5,3),fn_nPoints(5,2)) = -1;
coMatrix(fn_nPoints(5,3),fn_nPoints(5,3)) = 3+alpha;
coMatrix(fn_nPoints(5,3),fn_nPoints(4,3)) = -1;                          

%U(4,4)
coMatrix(fn_nPoints(4,4),fn_nPoints(4,4)) = 2+alpha*sqrt(2);
coMatrix(fn_nPoints(4,4),fn_nPoints(3,4)) = -1;
coMatrix(fn_nPoints(4,4),fn_nPoints(4,3)) = -1;

%U(3,5)
coMatrix(fn_nPoints(3,5),fn_nPoints(3,5)) = 3+alpha;
coMatrix(fn_nPoints(3,5),fn_nPoints(3,4)) = -1;
coMatrix(fn_nPoints(3,5),fn_nPoints(2,5)) = -1;

%% B Vector

%U(i,j) : i=0:2, j = 5

for i = 0:2
    coB(fn_nPoints(i,5)) = 50
end

%U(3,5)
coB(fn_nPoints(3,5)) = 50+alpha*tInf;

%U(4,4)
coB(fn_nPoints(4,4)) = alpha*sqrt(2)/2;

%U(5,3)
coB(fn_nPoints(5,3)) = alpha*tInf;
    
    
    