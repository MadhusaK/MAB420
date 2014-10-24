function genVisual(coMatrix,b,sol)

% A script to visualise the solution to the case study given in the lectures
% pretending that there is a piece missing, like in the assignment.

% Set up the parameters
n = 6;

% Generate the mesh
[X,Y] = meshgrid(0:1/n:1, 1:-1/n:0);

% Find the solution
% Tvec = coMatrix\b;
Tvec = sol;

% Load the solution into a matrix
T = zeros(size(X));

% Now pretend we have some bits missing
T(1,5:7) = nan;
T(2,5:7) = nan;
T(3,6:7) = nan;

% Load rows
T(1,1:4) = 0;
T(4:7, 7) = 50;

% load row 2
T(2,1:4) = Tvec(30:33);

% load row 3
T(3,1:5) = Tvec(25:29);


% load row 4
T(4,1:6) = Tvec(19:24);

% load row 5
T(5,1:6) = Tvec(13:18);


% load row 6
T(6,1:6) = Tvec(7:12);



% load row 7
T(7,1:6) = Tvec(1:6);



% Now extract the actual values we want
Tvec = T(~isnan(T));
Xvec = X(~isnan(T));
Yvec = Y(~isnan(T));

% Triangulate the domain
tri = delaunay(Xvec, Yvec);

% Check the triangulation
% figure
% trimesh(tri, Xvec, Yvec);
% axis image;

% Oops!  We have some extra triangles we don't want.
% Here's a saved triangulation with those triangles removed
load tri;

% Check the triangulation again
% figure
% trimesh(tri, Xvec, Yvec);
% axis image;

% Good, now we can visualise the solution
figure
trisurf(tri, Xvec, Yvec, Tvec);
title('Visualization of solution - Fullpacked Solution')
view(2);
shading interp;
colormap hot;
colorbar;
axis image;
grid off;
whitebg(gcf, [0.6 0.6 0.6]);
