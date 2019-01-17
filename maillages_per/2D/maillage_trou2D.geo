pas_mesh=0.01;
cx=0.5;
cy=0.5;
rayon=0.35;
xmin=0.;
ymin=0.;
xmax=1.;
ymax=1.;
Point(1) = {cx, cy, 0, pas_mesh};
Point(2) = {cx-rayon, cy, 0, pas_mesh};
Point(3) = {cx+rayon, cy, 0, pas_mesh};
Point(4) = {xmin, ymin, 0, pas_mesh};
Point(5) = {xmax, ymin, 0, pas_mesh};
Point(6) = {xmax, ymax, 0, pas_mesh};
Point(7) = {xmin, ymax, 0, pas_mesh};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 2};
// coté ouest
Line(3) = {7, 4};
// coté sud
Line(4) = {4, 5};
// coté est
Line(5) = {5, 6};
// coté nord
Line(6) = {6, 7};
// Forme le carré
Line Loop(7) = {3, 4, 5, 6};
// Forme le cercle
Line Loop(8) = {2, 1};
// Crée la surface à mailer (carree moins cercle)
Plane Surface(9) = {7, 8};
// Creation des lignes pour les conditions aux limites 
// 10 : cercle, 11 : ouest, 12: est, 13 :sud, 14 : nord
Physical Line(10) = {2, 1};
Physical Line(11) = {3};
Physical Line(12) = {5};
Physical Line(13) = {4};
Physical Line(14) = {6};
// On note 1 la surface complete
Physical Surface(1) = {9};
// definition de la periodicite face 3 et 5 (ouest/est)
Periodic Line {5}={-3};
// definition de la periodicite face 4 et 6 (sud/nord)
Periodic Line {6}={-4};

