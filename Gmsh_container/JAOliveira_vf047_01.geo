//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1, 1, 1};
//+
Radius = DefineNumber[ 0.3868, Name "Parameters/Radius" ];
//+
Cylinder(2) = {0, 0, 0, 1, 0, 0, Radius, 2*Pi};
//+
Cylinder(3) = {0, 1, 1, 1, 0, 0, Radius, 2*Pi};
//+
Cylinder(4) = {0, 0, 1, 1, 0, 0, Radius, 2*Pi};
//+
Cylinder(5) = {0, 1, 0, 1, 0, 0, Radius, 2*Pi};
//+
BooleanIntersection{ Volume{1}; }{ Volume{4}; Volume{3}; Volume{5}; Volume{2}; Delete; }
//+
Physical Volume("fiber", 1) = {4, 3, 5, 2};
//+
BooleanDifference{ Volume{1}; Delete; }{ Volume{4}; Volume{3}; Volume{5}; Volume{2}; }
//+
Physical Volume("matrix", 2) = {1};
//+
Physical Surface("top", 1) = {24, 30, 19};
//+
Physical Surface("bottom", 2) = {8, 28, 13};
//+
Physical Surface("left", 3) = {25, 31, 10};
//+
Physical Surface("right", 4) = {18, 29, 14};
//+
Physical Surface("front", 5) = {7, 12, 17, 22, 27};
//+
Physical Surface("back", 6) = {11, 16, 21, 26, 32};

//+
Split Curve {50} Point 	{};
//+
Split Curve {50} Point {};
//+
Point(25) = {0.0, 0.5, 1.0, 1.0};
//+
Point(26) = {1.0, 0.5, 1.0, 1.0};
//+
Point(27) = {1.0, 0.5, 0.0, 1.0};
//+
Point(28) = {0.0, 0.5, 0.0, 1.0};
//+
Line(53) = {26, 27};
//+
Line(54) = {27, 28};
//+
Line(55) = {28, 25};
//+
Line(56) = {25, 26};
