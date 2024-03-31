//+
SetFactory("OpenCASCADE");
//+
lc = 1e-2;
//+
Circle(1) = {0, 0, 0, 0.34};
Circle(2) = {0, 0, 0, 0.5};
Curve Loop(1) = {1};
Curve Loop(2) = {2};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
//+
Rectangle(3) = {-0.5, 0., 0, .5, .5};
Rectangle(4) = {0., 0., 0, .5, .5};
Rectangle(5) = {-0.5, -0.5, 0, .5, .5};
Rectangle(6) = {0., -0.5, 0, .5, .5};


//+
//Plane Surface(1) = {1};
//Physical Surface(1) = {1};
//e() = Extrude {0, 0, 1}{ Surface{14}; };

//Plane Surface(2) = {2};
//Physical Surface(2) = {2};

//+
BooleanFragments{ Surface{3}; Surface{4}; Surface{5}; Surface{6}; Delete; }{ Surface{1}; Surface{2}; Delete; }
//BooleanFragments{ Surface{10}; Surface{11}; Surface{12}; Surface{13}; Delete; }{ Surface{20}; Surface{21}; Delete; }
//BooleanDifference{Surface{2}; Delete;}{Surface{1}; Delete;}


//+
//Periodic Surface {2} = {1} Translate {1, 0, 0};
Periodic Curve {17} = {-6} Translate {1, 0, 0};
//Periodic Curve {23} = {13} Translate {1, 0, 0};
//Periodic Curve {5} = {13} Translate {0, 1, 0};
//Periodic Curve {15} = {22} Translate {0, 1, 0};
//
MeshSize {:} = 0.075;
//MeshSize {3,4,12,9} = 0.05; // fiber
//MeshSize {2,6,11,8} = 0.05; // inner matrix
//MeshSize {1,5,10,7} = 0.05; // outer matrix
//+
Physical Surface(1) = {3,4,9,12}; // fiber
Physical Surface(2) = {1,5,11,8}; // inner matrix
Physical Surface(3) = {2,6,7,10}; // outter matrix
//+

