lc = 1;

Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {0, 1, 0, lc};
Point(4) = {1, 1, 0, lc};
Point(5) = {0, 0, 1, lc};
Point(6) = {1, 0, 1, lc};
Point(7) = {0, 1, 1, lc};
Point(8) = {1, 1, 1, lc};

Line(1)  = {1,2};
Line(2)  = {2,4};
Line(3)  = {4,3};
Line(4)  = {3,1};
Line(5)  = {5,6};
Line(6)  = {6,8};
Line(7)  = {8,7};
Line(8)  = {7,5};
Line(9)  = {1,5};
Line(10) = {2,6};
Line(11) = {3,7};
Line(12) = {4,8};

Line Loop(1) = {1,10,-5,-9};
Line Loop(2) = {11,-7,-12,3};
Line Loop(3) = {12,-6,-10,2};
Line Loop(4) = {1,2,3,4};
Line Loop(5) = {4,9,-8,-11};
Line Loop(6) = {8,5,6,7};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};

Surface Loop(1) = {6, 5, 4, 1, 3, 2};
Volume(1) = {1};

Physical Surface("inflow") = {4};
Physical Surface("outflow") = {6};
Physical Surface("environment") = {5, 1, 3, 2};
Physical Volume("flow") = {1};

Transfinite Line "*" = 20 Using Bump 0.25;
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";