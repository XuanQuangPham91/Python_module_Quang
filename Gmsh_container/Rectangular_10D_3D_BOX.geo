//+
SetFactory("OpenCASCADE");
//+
d = 0.5;
r = 0.34;
lc = 1e-2;
// POINTS =====================================================================
Point(1) = {0, 0, 0, lc};

//=============================================================================
//=============================================================================
//Box(10) = {-d,-d,-d, 	2*d, 2*d, 2*d};
Box(2) = {-d, -r,-r,  	2*d, 2*r, 2*r};
Cylinder(1) = {-d, 0, 0, 2*d, 0, 0, r, 2*Pi};

Box(3) = {-d,-d,-d, 	2*d, (d-r), (d-r)};
Box(4) = {-d,-d,-r, 	2*d, d-r, 2*r};
Box(5) = {-d,-d,r, 	2*d, (d-r), (d-r)};
Box(6) = {-d,-r,r, 	2*d, 2*r, d-r};
Box(7) = {-d,r,r, 	2*d, (d-r), (d-r)};
Box(8) = {-d,r,-r, 	2*d, d-r, 2*r};
Box(9) = {-d,r,-d, 	2*d, (d-r), (d-r)};
Box(10) = {-d,-r,-d, 	2*d, 2*r, d-r};
//Box(6) = {r,r,-d, (d-r), (d-r), 2*d};
Box(50) = {-d, -d, -d, 2*d, 2*d, 2*d};
v() = BooleanFragments { Volume{50}; Delete; }{ Volume{1:10}; Delete; };
Geometry.OCCBoundsUseStl = 1;
eps = 1e-3; 
vin() = Volume In BoundingBox {-d-eps, -d-eps, -d-eps, 2*d+eps, 2*d+eps, 2*d+eps};
v() -= vin();


//Recursive Delete{ Volume{v()}; }
//+

//+Radius = DefineNumber[ 0.3868, Name "Parameters/Radius" ];
//+
//+Cylinder(9) = {0, 0, 0, 1, 0, 0, Radius, 2*Pi};
//+
//+Cylinder(10) = {0, 1, 1, 1, 0, 0, Radius, 2*Pi};
//+
//+Cylinder(11) = {0, 0, 1, 1, 0, 0, Radius, 2*Pi};
//+
//+Cylinder(12) = {0, 1, 0, 1, 0, 0, Radius, 2*Pi};
//+
//BooleanIntersection{ Volume{5}; Volume{1}; Volume{7}; Volume{3}; Volume{4}; Volume{8}; Volume{2}; Volume{6};  }{ Volume{10}; Volume{12}; Volume{9}; Volume{11}; Delete; }
//BooleanIntersection{ Volume{1};}{ Volume{2}; Delete; }
//+
//+BooleanDifference{ Volume{3}; Volume{7}; Volume{8}; Volume{4}; Volume{2}; Volume{6}; Volume{1}; Volume{5}; Delete; }{ Volume{10}; Volume{9}; Volume{11}; Volume{12}; Volume{13}; Volume{14}; Volume{15}; Volume{16}; }
//+
//BooleanIntersection{ Volume{1};}{ Volume{2}; Delete; }
//BooleanIntersection{ Volume{2};}{ Volume{1}; Delete;}
//BooleanDifference{ Volume{2}; Delete;}{ Volume{1}; }

//BooleanDifference{ Volume{3}; Volume{7}; Volume{8}; Volume{4}; Volume{2}; Volume{6}; Volume{1}; Volume{5}; Delete; }{ Volume{10}; Volume{9}; Volume{11}; Volume{12}; Volume{13}; Volume{14}; Volume{15}; Volume{16}; }
//+Physical Volume("fiber", 1) = {9, 10, 11, 12, 13, 14, 15, 16};

//=============================================================================
//Subdomain 1: x direction
Periodic Surface {49} = {50} Translate {1, 0, 0};
//=============================================================================
//Subdomain 2: x direction
Periodic Surface {15} = {11} Translate {1, 0, 0};
Periodic Surface {10} = {6} Translate {1, 0, 0};
Periodic Surface {20} = {16} Translate {1, 0, 0};
Periodic Surface {5} = {1} Translate {1, 0, 0};
//=============================================================================
//Subdomain 3 4 5 6 7 8 9 10: x direction
Periodic Surface {21} = {22} Translate {1, 0, 0};
Periodic Surface {27} = {28} Translate {1, 0, 0};
Periodic Surface {35} = {36} Translate {1, 0, 0};
Periodic Surface {46} = {45} Translate {1, 0, 0};
Periodic Surface {56} = {55} Translate {1, 0, 0};
Periodic Surface {52} = {51} Translate {1, 0, 0};
Periodic Surface {41} = {40} Translate {1, 0, 0};
Periodic Surface {32} = {31} Translate {1, 0, 0};
//=============================================================================
// y direction
Periodic Surface {57} = {37} Translate {0, 1, 0};
Periodic Surface {53} = {29} Translate {0, 1, 0};
Periodic Surface {42} = {23} Translate {0, 1, 0};
//=============================================================================
// z direction
Periodic Surface {58} = {43} Translate {0, 0, 1};
Periodic Surface {48} = {34} Translate {0, 0, 1};
Periodic Surface {39} = {25} Translate {0, 0, 1};
//=============================================================================
//+
Physical Volume("fiber", 1) = {1};
Physical Volume("matrix_2", 2) = {11, 12, 13, 14};
Physical Volume("matrix_3", 3) = {3};
Physical Volume("matrix_4", 4) = {4};
Physical Volume("matrix_5", 5) = {5};
Physical Volume("matrix_6", 6) = {6};
Physical Volume("matrix_7", 7) = {7};
Physical Volume("matrix_8", 8) = {8};
Physical Volume("matrix_9", 9) = {9};
Physical Volume("matrix_10", 10) = {10};



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
