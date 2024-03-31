//+

SetFactory("OpenCASCADE");
//
lc = 1e-2;
//

Point(1) = {-0.5, 0.5, 0, lc};
Point(2) = {0, 0.5, 0, lc};
Point(3) = {0.5, 0.5, 0, lc};
Point(4) = {-0.5, 0., 0, lc};
Point(5) = {0, 0., 0, lc};
Point(6) = {0.5, 0., 0, lc};
Point(7) = {-0.5, -0.5, 0, lc};
Point(8) = {0, -0.5, 0, lc};
Point(9) = {0.5, -0.5, 0, lc};

Rectangle(30) = {-2, -2, 0, 4, 4};

// point: 1 3 9 7
//Line(1) = {1, 3};
//Line(2) = {3, 9};
//Line(3) = {9, 7};
//Line(4) = {7, 1};
//Line Loop(1) = {1, 2, 3, 4};

//SetFactory("OpenCASCADE");

//-----------------------------------------------------------------------------
// point: 1 2 5 4
Line(1) = {1, 2};
Line(2) = {2, 5};
Line(3) = {5, 4};
Line(4) = {4, 1};
//Line Loop(1) = {1, 2, 5, 4};
// point: 2 3 6 5
Line(5) = {2, 3};
Line(6) = {3, 6};
Line(7) = {5, 6};
//Line Loop(2) = {5, 6, 7, 2};
// point: 5 6 9 8
Line(8) = {6, 9};
Line(9) = {9, 8};
Line(10) = {8, 5};
//Line Loop(3) = {7, 8, 9, 10};
// point: 4 5 8 7
Line(11) = {8, 7};
Line(12) = {7, 4};
//Line Loop(4) = {3, 10, 11, 12};
//Rectangle(4) = {.2, 0, .1, .3};
//+
Curve Loop(1) = {1, 5, 6, 8, 9, 11, 12, 4};
//Radius = DefineNumber[ 0.345, Name "Parameters/Radius" ];

//-----------------------------------------------------------------------------
//+ BooleanIntersection{ Volume{5}; Volume{1}; Volume{7}; Volume{3}; Volume{4}; Volume{8}; Volume{2}; Volume{6};  }{ Volume{10}; Volume{12}; Volume{9}; Volume{11}; Delete; }
//+
//+ BooleanDifference{ Volume{3}; Volume{7}; Volume{8}; Volume{4}; Volume{2}; Volume{6}; Volume{1}; Volume{5}; Delete; }{ Volume{10}; Volume{9}; Volume{11}; Volume{12}; Volume{13}; Volume{14}; Volume{15}; Volume{16}; }
//+
//+Physical Volume("fiber", 1) = {9, 10, 11, 12, 13, 14, 15, 16};

//+
//+Physical Volume("matrix", 2) = {1, 2, 3, 4, 5, 6, 7, 8};

//+
//+Physical Surface("top", 1) = {60, 65, 87, 92, 97, 102, 68, 72};
//+
//+Physical Surface("bottom", 2) = {52, 57, 116, 119, 112, 108, 79, 83};
//+
//+Physical Surface("left", 3) = {64, 59, 86, 91, 115, 118, 55, 50};
//+
//+Physical Surface("right", 4) = {69, 73, 98, 103, 107, 111, 78, 82};
//+
//+Physical Surface("front", 5) = {62, 75, 95, 100, 120, 113, 53, 84};
//+
//+Physical Surface("back", 6) = {63, 67, 85, 101, 105, 76, 114, 54};

//+
//Line(13) = {2, 5};
//Line(14) = {6, 4};
//Line(15) = {2, 8};
//Line(16) = {6, 4};
//+

//SetFactory("OpenCASCADE");
//
Circle(20) = {0, 0, 0, 0.345, 0, 2*Pi};
Curve Loop(2) = {20};
//
Circle(21) = {0, 0, 0, 0.48, 0, 2*Pi};
Curve Loop(3) = {21};
//+
//Plane Surface(1) = {1};
//Physical Surface(1) = {1};
e() = Extrude {0, 0, 0}{ Surface{10}; };

Plane Surface(2) = {2};
Physical Surface(2) = {2};

//+
//BooleanDifference{Surface{2}; Delete;}{Surface{1}; Delete;}


//+
//Periodic Surface {2} = {1} Translate {1, 0, 0};
Periodic Curve {6} = {4} Translate {1, 0, 0};
Periodic Curve {8} = {12} Translate {1, 0, 0};
Periodic Curve {1} = {5} Translate {0, 1, 0};
Periodic Curve {11} = {9} Translate {0, 1, 0};
//
MeshSize {:} = 0.1;
MeshSize {1} = 0.02;
//+

