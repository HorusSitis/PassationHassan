////////////// Pas du maillage //////////////

pas_mesh=0.01;

////////////// Paramètres géométriques //////////////

r_1=0.05;
r_2=0.1;
r_3=0.15;
r_4=0.2;
r_5=0.25;
r_6=0.3;
r_7=0.35;
r_8=0.4;
// extrapolation ?
r_9=0.45;

xmin=0.;
ymin=0.;
xmax=1.;
ymax=1.;

cx=0.5;
cy=0.5;

////////////// Indexation //////////////

npt_sq=0;
npt_cer1=10;
npt_cer2=20;
npt_cer3=30;
npt_cer4=40;
npt_cer5=50;
npt_cer6=60;
npt_cer7=70;
npt_cer8=80;

// lignes simples
nlext=100;
nlcer=200;
// boucles
nlo=1000;

// Surfaces planes
nSpl=3000;

// Lignes physiques
nPhl=4000;
// Surfaces physiques
nPhS=5000;

////////////// Construction du maillage fixe //////////////

/// Carré extérieur ///

Point(npt_sq+1) = {xmin, ymin, 0, pas_mesh};
Point(npt_sq+2) = {xmax, ymin, 0, pas_mesh};
Point(npt_sq+3) = {xmax, ymax, 0, pas_mesh};
Point(npt_sq+4) = {xmin, ymax, 0, pas_mesh};

// coté ouest
Line(nlext+1) = {npt_sq+4,npt_sq+1};
// coté sud
Line(nlext+2) = {npt_sq+1,npt_sq+2};
// coté est
Line(nlext+3) = {npt_sq+2,npt_sq+3};
// coté nord
Line(nlext+4) = {npt_sq+3,npt_sq+4};

// Forme le carré
Line Loop(nlo+nlext+1) = {nlext+1,nlext+2,nlext+3,nlext+4};

// definition de la periodicite face 1 et 3 (ouest/est)
Periodic Line {nlext+1}={-(nlext+3)};
// definition de la periodicite face 2 et 4 (sud/nord)
Periodic Line {nlext+2}={-(nlext+4)};

/// Cercles concentriques ///

//// r=0.05
Point(npt_cer1+1)={cx-r_1, cy, 0, pas_mesh};
Point(npt_cer1+2)={cx, cy, 0, pas_mesh};
Point(npt_cer1+3)={cx+r_1, cy, 0, pas_mesh};
Circle(nlcer+npt_cer1+1)={npt_cer1+3, npt_cer1+2, npt_cer1+1};
Circle(nlcer+npt_cer1+2)={npt_cer1+1, npt_cer1+2, npt_cer1+3};
Line Loop(nlo+nlcer+npt_cer1+1)={nlcer+npt_cer1+1, nlcer+npt_cer1+2};
//// r=0.1
Point(npt_cer2+1)={cx-r_2, cy, 0, pas_mesh};
Point(npt_cer2+2)={cx, cy, 0, pas_mesh};
Point(npt_cer2+3)={cx+r_2, cy, 0, pas_mesh};
Circle(nlcer+npt_cer2+1)={npt_cer2+3, npt_cer2+2, npt_cer2+1};
Circle(nlcer+npt_cer2+2)={npt_cer2+1, npt_cer2+2, npt_cer2+3};
Line Loop(nlo+nlcer+npt_cer2+1)={nlcer+npt_cer2+1, nlcer+npt_cer2+2};
//// r=0.15
Point(npt_cer3+1)={cx-r_3, cy, 0, pas_mesh};
Point(npt_cer3+2)={cx, cy, 0, pas_mesh};
Point(npt_cer3+3)={cx+r_3, cy, 0, pas_mesh};
Circle(nlcer+npt_cer3+1)={npt_cer3+3, npt_cer3+2, npt_cer3+1};
Circle(nlcer+npt_cer3+2)={npt_cer3+1, npt_cer3+2, npt_cer3+3};
Line Loop(nlo+nlcer+npt_cer3+1)={nlcer+npt_cer3+1, nlcer+npt_cer3+2};
//// r=0.2
Point(npt_cer4+1)={cx-r_4, cy, 0, pas_mesh};
Point(npt_cer4+2)={cx, cy, 0, pas_mesh};
Point(npt_cer4+3)={cx+r_4, cy, 0, pas_mesh};
Circle(nlcer+npt_cer4+1)={npt_cer4+3, npt_cer4+2, npt_cer4+1};
Circle(nlcer+npt_cer4+2)={npt_cer4+1, npt_cer4+2, npt_cer4+3};
Line Loop(nlo+nlcer+npt_cer4+1)={nlcer+npt_cer4+1, nlcer+npt_cer4+2};
//// r=0.25
Point(npt_cer5+1)={cx-r_5, cy, 0, pas_mesh};
Point(npt_cer5+2)={cx, cy, 0, pas_mesh};
Point(npt_cer5+3)={cx+r_5, cy, 0, pas_mesh};
Circle(nlcer+npt_cer5+1)={npt_cer5+3, npt_cer5+2, npt_cer5+1};
Circle(nlcer+npt_cer5+2)={npt_cer5+1, npt_cer5+2, npt_cer5+3};
Line Loop(nlo+nlcer+npt_cer5+1)={nlcer+npt_cer5+1, nlcer+npt_cer5+2};
//// r=0.3
Point(npt_cer6+1)={cx-r_6, cy, 0, pas_mesh};
Point(npt_cer6+2)={cx, cy, 0, pas_mesh};
Point(npt_cer6+3)={cx+r_6, cy, 0, pas_mesh};
Circle(nlcer+npt_cer6+1)={npt_cer6+3, npt_cer6+2, npt_cer6+1};
Circle(nlcer+npt_cer6+2)={npt_cer6+1, npt_cer6+2, npt_cer6+3};
Line Loop(nlo+nlcer+npt_cer6+1)={nlcer+npt_cer6+1, nlcer+npt_cer6+2};
//// r=0.35
Point(npt_cer7+1)={cx-r_7, cy, 0, pas_mesh};
Point(npt_cer7+2)={cx, cy, 0, pas_mesh};
Point(npt_cer7+3)={cx+r_7, cy, 0, pas_mesh};
Circle(nlcer+npt_cer7+1)={npt_cer7+3, npt_cer7+2, npt_cer7+1};
Circle(nlcer+npt_cer7+2)={npt_cer7+1, npt_cer7+2, npt_cer7+3};
Line Loop(nlo+nlcer+npt_cer7+1)={nlcer+npt_cer7+1, nlcer+npt_cer7+2};
//// r=0.4
Point(npt_cer8+1)={cx-r_8, cy, 0, pas_mesh};
Point(npt_cer8+2)={cx, cy, 0, pas_mesh};
Point(npt_cer8+3)={cx+r_8, cy, 0, pas_mesh};
Circle(nlcer+npt_cer8+1)={npt_cer8+3, npt_cer8+2, npt_cer8+1};
Circle(nlcer+npt_cer8+2)={npt_cer8+1, npt_cer8+2, npt_cer8+3};
Line Loop(nlo+nlcer+npt_cer8+1)={nlcer+npt_cer8+1, nlcer+npt_cer8+2};


/// Surfaces concentriques

Plane Surface(nSpl+nlo+nlcer+nptcer+1)={nlo+nlcer+npt_cer1+1};
Plane Surface(nSpl+nlo+nlcer+nptcer+2)={nlo+nlcer+npt_cer2+1, nlo+nlcer+npt_cer1+1};
Plane Surface(nSpl+nlo+nlcer+nptcer+3)={nlo+nlcer+npt_cer3+1, nlo+nlcer+npt_cer2+1};
Plane Surface(nSpl+nlo+nlcer+nptcer+4)={nlo+nlcer+npt_cer4+1, nlo+nlcer+npt_cer3+1};
Plane Surface(nSpl+nlo+nlcer+nptcer+5)={nlo+nlcer+npt_cer5+1, nlo+nlcer+npt_cer4+1};
Plane Surface(nSpl+nlo+nlcer+nptcer+6)={nlo+nlcer+npt_cer6+1, nlo+nlcer+npt_cer5+1};
Plane Surface(nSpl+nlo+nlcer+nptcer+7)={nlo+nlcer+npt_cer7+1, nlo+nlcer+npt_cer6+1};
Plane Surface(nSpl+nlo+nlcer+nptcer+8)={nlo+nlcer+npt_cer8+1, nlo+nlcer+npt_cer7+1};
//Plane Surface(nSpl+nlo+nlcer+nptcer+8)
Plane Surface(nSpl+nlo+nlcer+nptcer+10)={nlo+nlext+1, nlo+nlcer+npt_cer8+1};


// Creation des lignes  physiques : nécessaire ? Pas de condition aux limites.
// 11 : ouest, 12: est, 13 :sud, 14 : nord
//Physical Line(7) = {1};
//Physical Line(8) = {2};
//Physical Line(9) = {3};
//Physical Line(10) = {4};

// Surfaces physiques àmailler séparément
Physical Surface(nPhS+1) = {nSpl+nlo+nlcer+nptcer+1};
Physical Surface(nPhS+2) = {nSpl+nlo+nlcer+nptcer+2};
Physical Surface(nPhS+3) = {nSpl+nlo+nlcer+nptcer+3};
Physical Surface(nPhS+4) = {nSpl+nlo+nlcer+nptcer+4};
Physical Surface(nPhS+5) = {nSpl+nlo+nlcer+nptcer+5};
Physical Surface(nPhS+6) = {nSpl+nlo+nlcer+nptcer+6};
Physical Surface(nPhS+7) = {nSpl+nlo+nlcer+nptcer+7};
Physical Surface(nPhS+8) = {nSpl+nlo+nlcer+nptcer+8};
//Physical Surface(nPhS+9) = {nSpl+nlo+nlcer+nptce+9};
Physical Surface(nPhS+10) = {nSpl+nlo+nlcer+nptcer+10};


