// Mesh spacing
lc = 300;

// Points
Point(0) = {0.0,0.0,0,lc};
Point(1) = {60000.0,0.0,0,lc};
Point(2) = {60000.0,20000.0,0,lc};
Point(3) = {0.0,20000.0,0,lc};

// Lines
Line(0) = {0,1};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,0};

// Line loop
Line Loop(5) = {0,1,2,3};

// Surface
Plane Surface(6) = {5};

