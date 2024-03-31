//+
SetFactory("OpenCASCADE");
//+
lc = 1e-2;
//+

r=0.34;
d=0.5;

// POINTS =====================================================================
Point(1) = {0, 0, 0, lc};

Point(2) = {-d, -d, 0, lc};	Point(3) = {-d, -r, 0, lc};
Point(4) = {-r, -r, 0, lc};	Point(5) = {-r, -d, 0, lc};

Point(6) = {-d, +r, 0, lc};	Point(7) = {-d, +d, 0, lc};
Point(8) = {-r, +d, 0, lc};	Point(9) = {-r, +r, 0, lc};

Point(10) = {+r, +r, 0, lc};	Point(11) = {+r, +d, 0, lc};
Point(12) = {+d, +d, 0, lc};	Point(13) = {+d, +r, 0, lc};

Point(14) = {+r, -d, 0, lc};	Point(15) = {+r, -r, 0, lc};
Point(16) = {+d, -r, 0, lc};	Point(17) = {+d, -d, 0, lc};

Point(18) = {-r, +0, 0, lc};	Point(19) = {+0, +r, 0, lc};
Point(20) = {+r, +0, 0, lc};	Point(21) = {+0, -r, 0, lc};

// LINES ======================================================================
Line(2) = {2, 3};	Line(3) = {3, 4};
Line(4) = {4, 5};	Line(5) = {5, 2};

Line(6) = {6, 7};	Line(7) = {7, 8};
Line(8) = {8, 9};	Line(9) = {9, 6};

Line(10) = {10, 11};	Line(11) = {11, 12};
Line(12) = {12, 13};	Line(13) = {13, 10};

Line(14) = {14, 15};	Line(15) = {15, 16};
Line(16) = {16, 17};	Line(17) = {17, 14};

Line(18) = {4, 18};	Line(19) = {18, 9};
Line(20) = {9, 19};	Line(21) = {19, 10};
Line(22) = {10, 20};	Line(23) = {20, 15};
Line(24) = {15, 21};	Line(25) = {21, 4};

Line(26) = {3, 6};	Line(27) = {8, 11};
Line(28) = {13, 16};	Line(29) = {14, 5};


Circle(30) = {18,1,19};	Circle(31) = {19,1,20};
Circle(32) = {20,1,21};	Circle(33) = {21,1,18};


Curve Loop(35) = {2, 3, 4, 5};
Plane Surface(3) = {35};

Curve Loop(36) = {3, 26, 9, 19, 18};
Plane Surface(4) = {36};

Curve Loop(38) = {6, 7, 8, 9};
Plane Surface(5) = {38};

Curve Loop(39) = {8, 27, 10, 21, 20};
Plane Surface(6) = {39};

Curve Loop(41) = {10, 11, 12, 13};
Plane Surface(7) = {41};

Curve Loop(42) = {15, 23, 22, 13, 28};
Plane Surface(8) = {42};

Curve Loop(44) = {14, 15, 16, 17};
Plane Surface(9) = {44};

Curve Loop(45) = {4, 25, 24, 14, 29};
Plane Surface(10) = {45};

Curve Loop(34) = {30, 31, 32, 33};
Plane Surface(1) = {34};

Curve Loop(47) = {18, 19, 20, 21, 22, 23, 24, 25};
Plane Surface(11) = {47};
//Curve Loop(40) = {30, 31, 32, 33};
//Plane Surface(7) = {40};


//Circle(14) = {0, 0, 0, 0.382};
//Curve Loop(14) = 14;
//Plane Surface(20) = {14};
//+
//



//+
//Plane Surface(1) = {1};
Physical Surface(1) = {1};
Physical Surface(2) = {11,12,13,14};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
Physical Surface(6) = {6};
Physical Surface(7) = {7};
Physical Surface(8) = {8};
Physical Surface(9) = {9};
Physical Surface(10) = {10};

//e() = Extrude {0, 0, 1}{ Surface{14}; };

//Plane Surface(2) = {2};
//Physical Surface(2) = {2};

//+
//BooleanFragments{ Surface{10}; Surface{11}; Surface{12}; Surface{13}; Delete; }{ Surface{20}; Delete; }
//BooleanFragments{ Surface{10}; Surface{11}; Surface{12}; Surface{13}; Delete; }{ Surface{20}; Surface{21}; Delete; }
//BooleanDifference{Surface{2}; Delete;}{Surface{1}; Delete;}


//+
Periodic Surface {2} = {16} Translate {1, 0, 0};
Periodic Surface {26} = {28} Translate {1, 0, 0};
Periodic Surface {6} = {12} Translate {1, 0, 0};
//Periodic Curve {2} = {16} Translate {1, 0, 0};
//Periodic Curve {10} = {20} Translate {1, 0, 0};
//Periodic Curve {19} = {13} Translate {0, 1, 0};
//Periodic Curve {11} = {1} Translate {0, 1, 0};
//
MeshSize {:} = 0.025;

//Mesh.Algorithm = 6; // Frontal-Delaunay for 2D meshes
//MeshAlgorithm Surface {1} = 1; // MeshAdapt on surfaces 31 and 35

//MeshSize {1,11,12,13,14} = 0.05;
//MeshSize {2:10} = 0.05;
//MeshSize {3,4,12,9} = 0.05; // fiber
//MeshSize {2,6,11,8} = 0.05; // inner matrix
//MeshSize {1,5,10,7} = 0.05; // outer matrix
//+
//Physical Surface(1) = {7}; // fiber
//Physical Surface(2) = {8,9,10,11}; // inner matrix
//Physical Surface(3) = {1}; // matrix
//Physical Surface(4) = {2}; // matrix
//Physical Surface(5) = {2}; // matrix
//Physical Surface(6) = {2}; // matrix
//+

//+
Coherence;
