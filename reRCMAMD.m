function reRCMAMD(coMatrix,cholFull)
% 1. Computes the RCM and AMD node reordering for a given positive definite
% coefficient matrix, coMatrix
% 2. Computes Cholesky Factorization for both the RCM and AMD node
% reordering

RCM = symrcm(coMatrix);
AMD = symamd(coMatrix);

Arcm = coMatrix(RCM,RCM);
Aamd = coMatrix(AMD,AMD);

cholRCM = chol(Arcm);
cholAMD = chol(Aamd);

figure
subplot(2,3,1), spy(coMatrix), title('Coefficient Matrix');
subplot(2,3,2), spy(Arcm); title('RCM Node reodering');
subplot(2,3,3), spy(Aamd); title('AMD Node reodering');
subplot(2,3,4), spy(cholFull); title('Cholesky Full Packed');
subplot(2,3,5), spy(cholRCM); title('Cholesky RCM');
subplot(2,3,6), spy(cholAMD); title('Cholesky AMD');


