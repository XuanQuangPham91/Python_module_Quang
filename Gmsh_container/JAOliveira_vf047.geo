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
//+Physical Volume("fiber", 1) = {9, 10, 11, 12, 13, 14, 15, 16};

//+Physical Volume("matrix", 2) = {1, 2, 3, 4, 5, 6, 7, 8};

//+Physical Surface("top", 1) = {60, 65, 87, 92, 97, 102, 68, 72};
//+Physical Surface("bottom", 2) = {52, 57, 116, 119, 112, 108, 79, 83};
//+Physical Surface("left", 3) = {64, 59, 86, 91, 115, 118, 55, 50};
//+Physical Surface("right", 4) = {69, 73, 98, 103, 107, 111, 78, 82};
//+Physical Surface("front", 5) = {62, 75, 95, 100, 120, 113, 53, 84};
//+Physical Surface("back", 6) = {63, 67, 85, 101, 105, 76, 114, 54};



//+
Line(210) = {101, 63};
//+
Line(109) = {38, 39};
//+
Line(110) = {45, 63};
//+
Line(111) = {27, 3};
//+
Line(112) = {3, 2};
//+
Line(113) = {48, 45};
//+
Line(114) = {12, 9};
//+
Line(115) = {17, 3};
//+
Line(116) = {53, 45};
//+
Line(203) = {38, 101};
//+
Line(204) = {95, 2};
//+
Line(205) = {17, 95};
//+
Line(206) = {101, 53};
//+
Line(207) = {27, 95};
//+
Line(208) = {101, 48};
//+
Line(209) = {12, 95};
//+
Curve Loop(121) = {100, 203, -172, -165, 204, 109};
//+
Plane Surface(121) = {121};
//+
Curve Loop(122) = {148, 154, -208, -172, -165, -209};
//+
Plane Surface(122) = {122};
//+
Curve Loop(124) = {125, 117, -206, -172, -165, -205};
//+
Plane Surface(123) = {124};
//+
Curve Loop(126) = {165, 172, 210, -136, -130, 207};
//+
Plane Surface(124) = {126};
//+ BooleanFragments{ Volume{1}; Volume{10}; Volume{5}; Volume{9}; Volume{15}; Volume{2}; Volume{6}; Volume{16}; Volume{8}; Volume{14}; Volume{4}; Volume{13}; Volume{12}; Volume{3}; Volume{11}; Volume{7}; Delete;}{ Surface{124}; Surface{123}; Surface{122}; Surface{121}; Delete; }

//+BooleanUnion{ Volume{28}; Volume{27}; Volume{32}; Volume{31}; Volume{22}; Volume{18}; Volume{21}; Volume{17}; Delete; }{ Volume{1}; Volume{5}; Volume{2}; Volume{6}; Volume{11}; Volume{13}; Volume{14}; Volume{12}; Delete; }

//+BooleanUnion{ Volume{26}; Volume{30}; Volume{23}; Volume{19}; Volume{16}; Volume{10}; Volume{7}; Volume{3}; Delete; }{ Volume{25}; Volume{29}; Volume{24}; Volume{20}; Volume{9}; Volume{15}; Volume{4}; Volume{8}; Delete; }
//+ BooleanUnion{ Volume{25}; Volume{26}; Volume{24}; Volume{23}; Volume{9}; Volume{10}; Volume{3}; Volume{4}; Delete; }{ Volume{19}; Volume{20}; Volume{15}; Volume{16}; Volume{30}; Volume{29}; Volume{7}; Volume{8}; Delete; }
