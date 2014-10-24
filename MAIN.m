% Main Script
close all
[coMatrix,b] = genAB(); % Generate coefficient and b matrix

% Uncomment whichever code you wish to run

[solFull,cholFac,cholFac2] = storFull(coMatrix,b);    % Full Packed Solution
genVisual(coMatrix,b,solFull)                         % Visualisation
reRCMAMD(coMatrix,cholFac);                           % Node Reordering
[pacVector,cholPac,solPac] = storPac(coMatrix,b);     % Packed Solution
[banMatrix,cholBan,solBan] = storBan(coMatrix,b);     % Band Solution
[sparMatrix,cholSpar,solSpar] = storSpar(coMatrix,b); % Sparse Solution

