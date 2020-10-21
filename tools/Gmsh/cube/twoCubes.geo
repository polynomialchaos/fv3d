lc = 1;

Point(1) = {0, 0, -1, lc};
Point(2) = {1, 0, -1, lc};
Point(3) = {0, 1, -1, lc};
Point(4) = {1, 1, -1, lc};
Point(5) = {0, 0, 0, lc};
Point(6) = {1, 0, 0, lc};
Point(7) = {0, 1, 0, lc};
Point(8) = {1, 1, 0, lc};
Point(9) = {0, 0, 1, lc};
Point(10) = {1, 0, 1, lc};
Point(11) = {0, 1, 1, lc};
Point(12) = {1, 1, 1, lc};//+
Line(1) = {11, 7};
//+
Line(2) = {7, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 8};
//+
Line(5) = {8, 7};
//+
Line(6) = {8, 12};
//+
Line(7) = {12, 11};
//+
Line(8) = {9, 5};
//+
Line(9) = {5, 1};
//+
Line(10) = {1, 2};
//+
Line(11) = {2, 6};
//+
Line(12) = {6, 10};
//+
Line(13) = {10, 9};
//+
Line(14) = {5, 6};
//+
Line(15) = {9, 11};
//+
Line(16) = {5, 7};
//+
Line(17) = {1, 3};
//+
Line(18) = {12, 10};
//+
Line(19) = {6, 8};
//+
Line(20) = {2, 4};
//+
Line Loop(1) = {16, -5, -19, -14};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {17, 3, -20, -10};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {15, -7, 18, 13};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {1, -5, 6, 7};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {2, 3, 4, 5};
//+
Plane Surface(5) = {5};
//+
Line Loop(6) = {14, -11, -10, -9};
//+
Plane Surface(6) = {6};
//+
Line Loop(7) = {8, 14, 12, 13};
//+
Plane Surface(7) = {7};
//+
Line Loop(8) = {18, -12, 19, 6};
//+
Plane Surface(8) = {8};
//+
Line Loop(9) = {11, 19, -4, -20};
//+
Plane Surface(9) = {9};
//+
Line Loop(10) = {17, -2, -16, 9};
//+
Plane Surface(10) = {10};
//+
Line Loop(11) = {1, -16, -8, 15};
//+
Plane Surface(11) = {11};
//+
Physical Surface("inflow") = {2};
//+
Physical Surface("outflow") = {3};
//+
Physical Surface("environment") = {10, 6, 5, 9, 8, 4, 11, 7};
//+
Surface Loop(1) = {5, 10, 2, 9, 6, 4, 11, 7, 8, 3};
//+
//Volume(1) = {1};
//+
//Physical Volume("flow") = {1};


Surface Loop(2) = {10, 2, 5, 9, 6, 1};
//+
Volume(1) = {2};
//+
Surface Loop(3) = {1, 8, 3, 11, 4, 7};
//+
Volume(2) = {3};
//+
Physical Volume("flow") = {1, 2};

Transfinite Line "*" = 2 Using Bump 1;
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";