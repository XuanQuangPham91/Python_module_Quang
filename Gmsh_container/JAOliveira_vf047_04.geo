//+
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 0.5, 0.5, 0.5};
//+
Box(2) = {0, 0.5, 0, 0.5, 0.5, 0.5};
//+
Box(3) = {0, 0, 0.5, 0.5, 0.5, 0.5};
//+
Box(4) = {0, 0.5, 0.5, 0.5, 0.5, 0.5};
//+
Box(5) = {0.5, 0, 0, 0.5, 0.5, 0.5};
//+
Box(6) = {0.5, 0.5, 0, 0.5, 0.5, 0.5};
//+
Box(7) = {0.5, 0, 0.5, 0.5, 0.5, 0.5};
//+
Box(8) = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
//+

Radius = DefineNumber[ 0.3868, Name "Parameters/Radius" ];
//+
Cylinder(9) = {0, 0, 0, 1, 0, 0, Radius, 2*Pi};
//+
Cylinder(10) = {0, 1, 1, 1, 0, 0, Radius, 2*Pi};
//+
Cylinder(11) = {0, 0, 1, 1, 0, 0, Radius, 2*Pi};
//+
Cylinder(12) = {0, 1, 0, 1, 0, 0, Radius, 2*Pi};
//+
BooleanIntersection{ Volume{5}; Volume{1}; Volume{7}; Volume{3}; Volume{4}; Volume{8}; Volume{2}; Volume{6};  }{ Volume{10}; Volume{12}; Volume{9}; Volume{11}; Delete; }
//+
BooleanDifference{ Volume{3}; Volume{7}; Volume{8}; Volume{4}; Volume{2}; Volume{6}; Volume{1}; Volume{5}; Delete; }{ Volume{10}; Volume{9}; Volume{11}; Volume{12}; Volume{13}; Volume{14}; Volume{15}; Volume{16}; }
//+
Physical Volume("fiber", 1) = {9, 10, 11, 12, 13, 14, 15, 16};

//+
Physical Volume("matrix", 2) = {1, 2, 3, 4, 5, 6, 7, 8};

//+
Physical Surface("top", 1) = {60, 65, 87, 92, 97, 102, 68, 72};
//+
Physical Surface("bottom", 2) = {52, 57, 116, 119, 112, 108, 79, 83};
//+
Physical Surface("left", 3) = {64, 59, 86, 91, 115, 118, 55, 50};
//+
Physical Surface("right", 4) = {69, 73, 98, 103, 107, 111, 78, 82};
//+
Physical Surface("front", 5) = {62, 75, 95, 100, 120, 113, 53, 84};
//+
Physical Surface("back", 6) = {63, 67, 85, 101, 105, 76, 114, 54};




