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

subplot(2,3,1), spy(coMatrix);
subplot(2,3,2), spy(Arcm);
subplot(2,3,3), spy(Aamd);
subplot(2,3,4), spy(cholFull);
subplot(2,3,5), spy(cholRCM);
subplot(2,3,6), spy(cholAMD);


