//+
SetFactory("OpenCASCADE");
//+
lc = 1e-2;

// 4 corner and center
Point(1) = {0, 0, 0, lc};
Point(2) = {-0.5, -0.5, 0, lc};
Point(3) = {-0.5, 0.5, 0, lc}; 
Point(4) = {0.5, 0.5, 0, lc}; 
Point(5) = {0.5, -0.5, 0, lc}; 


// 4 midles 
Point(6) = {-0.5, 0, 0, lc};
Point(7) = {0, 0.5, 0, lc}; 
Point(8) = {0.5, 0, 0, lc}; 
Point(9) = {0, -0.5, 0, lc};

//+
Circle(1) = {-0.5, -0.5, 0, 0.382, Pi/2}; 
Circle(2) = {-0.5, 0.5, 0, 0.382, -Pi/2, 0}; 
Circle(3) = {0.5, 0.5, 0, 0.382, Pi, -Pi/2};  
Circle(4) = {0.5, -0.5, 0, 0.382, Pi/2, Pi};

// 4corners surfaces
Line(5) = {11, 2};
Line(6) = {2, 10};
Curve Loop(1) = {1,5,6};  Plane Surface(1) = {1};   
Line(7) = {3, 12};
Line(8) = {13, 3};
Curve Loop(2) = {2,8,7};  Plane Surface(2) = {2};
Line(9) = {14, 4};
Line(10) = {4, 15};
Curve Loop(3) = {3,9,10};  Plane Surface(3) = {3};
Line(11) = {17, 5};
Line(12) = {5, 16};
Curve Loop(5) = {4,11,12};  Plane Surface(5) = {5};
//
Physical Surface(1) = {1,2,3,5};
//+
Line(13) = {10, 9};
Line(14) = {9, 1};
Line(15) = {1, 6};
Line(16) = {6, 11};
Curve Loop(6) = {1,13,14,15,16};  Plane Surface(6) = {6};
Line(17) = {13, 7};
Line(18) = {7, 1};
Line(19) = {6, 12};
Curve Loop(8) = {2,17,18,15,19};  Plane Surface(8) = {8};
Line(20) = {15, 8};
Line(21) = {8, 1};
Line(22) = {7, 14};
Curve Loop(9) = {3,20,21,-18,22};  Plane Surface(9) = {9};
Line(23) = {17, 9};
Line(24) = {8, 16};
Curve Loop(10) = {4,23,14,-21,24};  Plane Surface(10) = {10};
//
Physical Surface(2) = {6,8,9,10};
//+
//Rectangle(10) = {-0.5, 0., 0, .5, .5};
//Rectangle(11) = {0., 0., 0, .5, .5};
//Rectangle(12) = {-0.5, -0.5, 0, .5, .5};
//Rectangle(13) = {0., -0.5, 0, .5, .5};


//+
//Plane Surface(1) = {1};
//Physical Surface(1) = {1};
//e() = Extrude {0, 0, 1}{ Surface{14}; };

//Plane Surface(2) = {2};
//Physical Surface(2) = {2};

//+
//BooleanFragments{ Surface{10}; Surface{11}; Surface{12}; Surface{13}; Delete; }{ Surface{1}; Surface{2}; Surface{3}; Surface{4}; Delete; }
//BooleanFragments{ Surface{10}; Surface{11}; Surface{12}; Surface{13}; Delete; }{ Surface{20}; Surface{21}; Delete; }
//BooleanDifference{Surface{2}; Delete;}{Surface{1}; Delete;}
//Delete {Surface{9,10,11,12};}
//Delete {Curve {25,26,27,28,29,30};}
//Delete {Point {18,19};}
//+
Periodic Surface {8} = {6} Translate {1, 0, 0};
Periodic Curve {17} = {13} Translate {1, 0, 0};
Periodic Curve {22} = {23} Translate {1, 0, 0};
Periodic Curve {9} = {11} Translate {1, 0, 0};


Periodic Curve {10} = {7} Translate {0, 1, 0};
Periodic Curve {20} = {19} Translate {0, 1, 0};
Periodic Curve {24} = {16} Translate {0, 1, 0};
Periodic Curve {12} = {5} Translate {0, 1, 0};

//
MeshSize {:} = 0.015;
//MeshSize {3,4,12,9} = 0.05; // fiber
//MeshSize {2,6,11,8} = 0.05; // inner matrix
//MeshSize {1,5,10,7} = 0.05; // outer matrix
//+
//Physical Surface(1) = {2,4,5,8}; // fiber
//Physical Surface(2) = {1,3,6,7}; // inner matrix
//+

