// Gmsh project created on Fri Jan 28 16:33:55 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 30, 0, 1.0};
//+
Point(4) = {0, 30, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Transfinite Curve {3, 1} = 9 Using Progression 1;
//+
Transfinite Curve {4, 2} = 61 Using Progression 1;
//+//+
Transfinite Surface {1};
