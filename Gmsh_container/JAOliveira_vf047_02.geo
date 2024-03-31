//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1, 0.5, 0.5};
//+
Box(2) = {0, 0.5, 0, 1, 0.5, 0.5};
//+
Box(3) = {0, 0, 0, 1, 0.5, 1};
//+
Box(4) = {0, 0.5, 0, 1, 0.5, 1};
//+
Radius = DefineNumber[ 0.3868, Name "Parameters/Radius" ];
//+
Cylinder(6) = {0, 0, 0, 1, 0, 0, Radius, 2*Pi};
//+
Cylinder(7) = {0, 1, 1, 1, 0, 0, Radius, 2*Pi};
//+
Cylinder(8) = {0, 0, 1, 1, 0, 0, Radius, 2*Pi};
//+
Cylinder(5) = {0, 1, 0, 1, 0, 0, Radius, 2*Pi};
//+
BooleanIntersection{ Volume{1}; Volume{2}; Volume{3}; Volume{4};}{ Volume{7}; Volume{8}; Volume{5}; Volume{6}; Delete; }
//+
BooleanDifference{ Volume{1}; Volume{2}; Volume{3}; Volume{4}; Delete;}{ Volume{7}; Volume{8}; Volume{5}; Volume{6}; }


//+
Physical Volume("fiber", 1) = {5, 6, 7, 8};

//+
Physical Volume("matrix", 2) = {9, 10, 2, 1};

//+
Physical Surface("top", 1) = {37, 47, 52, 41};
//+
Physical Surface("bottom", 2) = {28, 59, 63, 33};
//+
Physical Surface("left", 3) = {36, 46, 57, 26};
//+
Physical Surface("right", 4) = {42, 53, 62, 32};
//+
Physical Surface("front", 5) = {39, 44, 29, 34, 60, 64, 50, 55};
//+
Physical Surface("back", 6) = {35, 40, 45, 51, 56, 61, 25, 30};

//+
Transfinite Curve {57, 51, 60, 66, 65, 59, 49, 53, 71, 67, 73, 68, 75, 76, 78, 81} = 8 Using Progression 1;
//+
Transfinite Curve {90, 94, 97, 84, 88, 104, 101, 83, 99,  95,103,107, 108, 109, 110, 111} = 3 Using Progression 1;
//+
Transfinite Curve {92, 85, 86, 96, 100, 93, 106, 102} = 10 Using Progression 1;
//+
Transfinite Curve {48, 34, 64, 52, 62, 105, 55, 63, 98, 80, 72, 89, 79, 70, 87, 54, 91} = 19 Using Progression 1;
//+
Transfinite Curve {56, 74, 69, 50, 82, 77, 61, 58} = 12 Using Progression 1;
