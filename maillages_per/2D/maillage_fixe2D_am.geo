pas_mesh=0.01;

xmin=0.;
ymin=0.;
xmax=1.;
ymax=1.;

Point(1) = {xmin, ymin, 0, pas_mesh};
Point(2) = {xmax, ymin, 0, pas_mesh};
Point(3) = {xmax, ymax, 0, pas_mesh};
Point(4) = {xmin, ymax, 0, pas_mesh};

// coté ouest
Line(1) = {4,1};
// coté sud
Line(2) = {1,2};
// coté est
Line(3) = {2,3};
// coté nord
Line(4) = {3,4};

// Forme le carré
Line Loop(5) = {1,2,3,4};

// Crée la surface à mailler (carree)
Plane Surface(6) = {5};

// Creation des lignes pour les conditions aux limites 
// 11 : ouest, 12: est, 13 :sud, 14 : nord
Physical Line(7) = {1};
Physical Line(8) = {2};
Physical Line(9) = {3};
Physical Line(10) = {4};

// On note 1 la surface complete
Physical Surface(1) = {6};

// definition de la periodicite face 1 et 3 (ouest/est)
Periodic Line {1}={-3};
// definition de la periodicite face 2 et 4 (sud/nord)
Periodic Line {2}={-4};

