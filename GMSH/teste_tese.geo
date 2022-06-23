// Gmsh project created on Fri Oct 29 11:37:18 2021
c1_1 = 0.125;
Point(1) = {0, 0, 0, c1_1};
Point(2) = {1, 0, 0, c1_1};
Point(3) = {1, 1, 0, c1_1};
Point(4) = {0, 1, 0, c1_1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {4, 1, 2, 3};
Plane Surface(6) = {5};
Physical Point(7) = {1, 2, 3, 4};
Physical Line(8) = {1, 2, 3, 4};
Physical Surface(9) = {6};
