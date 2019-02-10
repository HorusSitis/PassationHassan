pas_mesh=0.01;

rayon=0.35;

xmin=0.;
ymin=0.;
xmax=1.;
ymax=1.;

// sommets du carré, contenus dans l'inclusion solide

Point(1) = {xmin, ymin, 0, pas_mesh};
Point(2) = {xmax, ymin, 0, pas_mesh};
Point(3) = {xmax, ymax, 0, pas_mesh};
Point(4) = {xmin, ymax, 0, pas_mesh};

// extrémités des arcs

Point(5)= {xmin, ymin+rayon, 0, pas_mesh};
Point(6)= {xmin+rayon, ymin, 0, pas_mesh};
Point(7)= {xmax-rayon, ymin, 0, pas_mesh};
Point(8)= {xmax, ymin+rayon, 0, pas_mesh};
Point(9)= {xmax, ymax-rayon, 0, pas_mesh};
Point(10)= {xmax-rayon, ymax, 0, pas_mesh};
Point(11)= {xmin+rayon, ymax, 0, pas_mesh};
Point(12)= {xmin, ymax-rayon, 0, pas_mesh};

// quarts de cercles aux sommets

Circle(1) = {5,1,6};
Circle(2) = {7,2,8};
Circle(3) = {9,3,10};
Circle(4) = {11,4,12};

// coté sud
Line(5) = {6,7};
// coté est
Line(6) = {8,9};
// coté nord
Line(7) = {10,11};
// coté ouest
Line(8) = {12,5};

// Forme la piece de puzzle
Line Loop(9) = {1,5,2,6,3,7,4,8};

// Crée la surface à mailler (carree moins cercle)
Plane Surface(10) = {9};

// Creation des lignes pour les conditions aux limites 
//// 11 : solide ; 12 : sud , 13 : est, 14 : nord, 15 : ouest
Physical Line(11) = {1,2,3,4};
Physical Line(12) = {5};
Physical Line(13) = {6};
Physical Line(14) = {7};
Physical Line(15) = {8};

// On note 1 la surface complete
Physical Surface(1) = {10};

// definition de la periodicite face 8 et 6 (ouest/est)
Periodic Line {8}={-6};
// definition de la periodicite face 5 et 7 (sud/nord)
Periodic Line {5}={-7};

