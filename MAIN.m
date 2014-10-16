% Main Script

[coMatrix,b] = genAB();

% testMatrix= full(gallery('poisson',20));
% b = zeros(20,1);

[solFull,cholFac,cholFac2] = storFull(coMatrix,b);
% [pacVector,cholPac,solPac] = storPac(coMatrix,b);
% [banMatrix,cholBan,solBan] = storBan(coMatrix,b);
% [sparMatrix,cholSpar] = storSpar(coMatrix);

