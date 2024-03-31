//+
SetFactory("OpenCASCADE");
//+
lc = 1e-2;
//+
Circle(14) = {0, 0, 0, 0.382};
Curve Loop(14) = 14;
Plane Surface(20) = {14};
//+
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
BooleanFragments{ Surface{10}; Surface{11}; Surface{12}; Surface{13}; Delete; }{ Surface{20}; Delete; }
//BooleanFragments{ Surface{10}; Surface{11}; Surface{12}; Surface{13}; Delete; }{ Surface{20}; Surface{21}; Delete; }
//BooleanDifference{Surface{2}; Delete;}{Surface{1}; Delete;}


//+
//Periodic Surface {2} = {1} Translate {1, 0, 0};
Periodic Curve {2} = {16} Translate {1, 0, 0};
Periodic Curve {10} = {20} Translate {1, 0, 0};
Periodic Curve {19} = {13} Translate {0, 1, 0};
Periodic Curve {11} = {1} Translate {0, 1, 0};
//
MeshSize {:} = 0.015;
//MeshSize {3,4,12,9} = 0.05; // fiber
//MeshSize {2,6,11,8} = 0.05; // inner matrix
//MeshSize {1,5,10,7} = 0.05; // outer matrix
//+
Physical Surface(1) = {2,3,6,8}; // fiber
Physical Surface(2) = {1,4,7,5}; // inner matrix
//+

//+
Coherence;
