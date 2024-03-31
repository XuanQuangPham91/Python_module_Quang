//+
SetFactory("OpenCASCADE");
//+
lc = 1e-2;
//+
Circle(1) = {-0.5, -0.5, 0, 0.382}; Curve Loop(1) = 1;
Plane Surface(1) = {1};
//+
Circle(2) = {-0.5, 0.5, 0, 0.382}; Curve Loop(2) = 2;
Plane Surface(2) = {2};
//+
Circle(3) = {0.5, 0.5, 0, 0.382};  Curve Loop(3) = 3;
Plane Surface(3) = {3};
//+
Circle(4) = {0.5, -0.5, 0, 0.382};
Curve Loop(4) = 4;
Plane Surface(4) = {4};
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
BooleanFragments{ Surface{10}; Surface{11}; Surface{12}; Surface{13}; Delete; }{ Surface{1}; Surface{2}; Surface{3}; Surface{4}; Delete; }
//BooleanFragments{ Surface{10}; Surface{11}; Surface{12}; Surface{13}; Delete; }{ Surface{20}; Surface{21}; Delete; }
//BooleanDifference{Surface{2}; Delete;}{Surface{1}; Delete;}
Delete {Surface{9,10,11,12};}
Delete {Curve {25,26,27,28,29,30};}
Delete {Point {18,19};}
//+
//Periodic Surface {2} = {1} Translate {1, 0, 0};
Periodic Curve {7} = {16} Translate {1, 0, 0};
Periodic Curve {3} = {19} Translate {1, 0, 0};
Periodic Curve {8} = {22} Translate {1, 0, 0};
Periodic Curve {12} = {24} Translate {1, 0, 0};

Periodic Curve {13} = {6} Translate {0, 1, 0};
Periodic Curve {10} = {1} Translate {0, 1, 0};
Periodic Curve {20} = {17} Translate {0, 1, 0};
Periodic Curve {23} = {14} Translate {0, 1, 0};
//
MeshSize {:} = 0.02;
//MeshSize {3,4,12,9} = 0.05; // fiber
//MeshSize {2,6,11,8} = 0.05; // inner matrix
//MeshSize {1,5,10,7} = 0.05; // outer matrix
//+
Physical Surface(1) = {2,4,5,8}; // fiber
Physical Surface(2) = {1,3,6,7}; // inner matrix
//+

