//+
SetFactory("OpenCASCADE");
//+
lc = 1e-2;
//+
Circle(14) = {0, 0, 0, 0.345};
Curve Loop(14) = 14;
Plane Surface(20) = {14};
//+
Circle(16) = {0, 0, 0, 0.5};
Curve Loop(17) = 16;
Plane Surface(21) = {17};

//+
//Point(1) = {-0.5, 0.5, 0, lc};
//Point(2) = {0, 0.5, 0, lc};
//Point(3) = {0.5, 0.5, 0, lc};
//Point(4) = {-0.5, 0., 0, lc};
//Point(5) = {0, 0., 0, lc};
//Point(6) = {0.5, 0., 0, lc};
//Point(7) = {-0.5, -0.5, 0, lc};
//Point(8) = {0, -0.5, 0, lc};
//Point(9) = {0.5, -0.5, 0, lc};

//+
//Rectangle(s(0)) = {xmin, ymin, zmin, L, H};
//Rectangle(30) = {-0.5, -0.5, 0, 1, 1};
Rectangle(10) = {-0.5, 0., 0, .5, .5};
Rectangle(11) = {0., 0., 0, .5, .5};
Rectangle(12) = {-0.5, -0.5, 0, .5, .5};
Rectangle(13) = {0., -0.5, 0, .5, .5};


//+
//Plane Surface(1) = {1};
//Physical Surface(1) = {1};
//e() = Extrude {0, 0, 1}{ Surface{14}; };

//Plane Surface(2) = {2};
//Physical Surface(2) = {2};

//+
BooleanFragments{ Surface{10}; Surface{11}; Surface{12}; Surface{13}; Delete; }{ Surface{20}; Surface{21}; Delete; }

//BooleanDifference{Surface{2}; Delete;}{Surface{1}; Delete;}


//+
//Periodic Surface {2} = {1} Translate {1, 0, 0};
Periodic Curve {6} = {17} Translate {1, 0, 0};
Periodic Curve {13} = {23} Translate {1, 0, 0};
Periodic Curve {14} = {5} Translate {0, 1, 0};
Periodic Curve {22} = {15} Translate {0, 1, 0};
//
MeshSize {:} = 0.04;
//MeshSize {3,4,12,9} = 0.05; // fiber
//MeshSize {2,6,11,8} = 0.05; // inner matrix
//MeshSize {1,5,10,7} = 0.05; // outer matrix
//+
Physical Surface(1) = {3,4,12,9}; // fiber
Physical Surface(2) = {1,5,11,8}; // inner matrix
Physical Surface(3) = {6,10,7,2}; // outer matrix
//+

