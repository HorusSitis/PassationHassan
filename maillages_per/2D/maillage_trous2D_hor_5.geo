pas_mesh=0.01;

rayon_centre=0.05;
rayon_lat=0.15;

xmin=0.;
ymin=0.;
xmax=1.;
ymax=1.;
xcen=0.5;
ycen=0.5;


// Numero (moins 1) du premier point
npt=100;
// Numero (moins 1) du premier arc de cercle
narc=200;
// Numero (moins 1) de la premiere ligne droite
nldr=300;
// Numero (moins 1) de la premiere ligne composée
nlco=400;
// Numero (moins 1) de la premiere ligne physique
nlPh=500;
// Numero (moins 1) de la surface à mailler
nsurf=600;
// Numéro (moins 1) de la premiere surface physique
nsPh=700;
// Numéro (moins 1) de la premiere ligne periodique
nlPer=800;

// sommets du domaine fluide

Point(npt+1) = {xmin, ymin, 0, pas_mesh};
Point(npt+2) = {xmax, ymin, 0, pas_mesh};
Point(npt+3) = {xmax, ycen-rayon_lat, 0, pas_mesh};
Point(npt+4) = {xmax, ycen+rayon_lat, 0, pas_mesh};
Point(npt+5)= {xmax, ymax, 0, pas_mesh};
Point(npt+6)= {xmin, ymax, 0, pas_mesh};
Point(npt+7)= {xmin, ycen+rayon_lat, 0, pas_mesh};
Point(npt+8)= {xmin, ycen-rayon_lat, 0, pas_mesh};

// points pour générer les arcs de cercle

Point(npt+9)= {xmin, ycen, 0, pas_mesh};
Point(npt+10)= {xmax, ycen, 0, pas_mesh};
Point(npt+11)= {xcen, ycen, 0, pas_mesh};
Point(npt+12)= {xcen+rayon_centre, ycen, 0, pas_mesh};
Point(npt+13)= {xcen-rayon_centre, ycen, 0, pas_mesh};

Point(npt+14)={xmax-rayon_lat,ycen,0,pas_mesh};
Point(npt+15)={xmin+rayon_lat,ycen,0,pas_mesh};

// demis cercle pour les deux inclusions solides

Circle(narc+1) = {npt+3,npt+10,npt+14};
Circle(narc+2) = {npt+14,npt+10,npt+4};
Circle(narc+3) = {npt+7,npt+9,npt+15};
Circle(narc+4) = {npt+15,npt+9,npt+8};

Circle(narc+5) = {npt+12,npt+11,npt+13};
Circle(narc+6) = {npt+13,npt+11,npt+12};

// cote sud
Line(nldr+1) = {npt+1,npt+2};
// cote nord
Line(nldr+4) = {npt+5,npt+6};

// cote est
Line(nldr+2)={npt+2,npt+3};
Line(nldr+3)={npt+4,npt+5};
Line Loop(nlco+1)={nldr+2,narc+1,narc+2,nldr+3};

// cote ouest
Line(nldr+5)={npt+6,npt+7};
Line(nldr+6)={npt+8,npt+1};
Line Loop(nlco+2)={nldr+5,narc+3,narc+4,nldr+6};

// Forme la piece de puzzle
Line Loop(nlco+3) = {nldr+1,nldr+2,narc+1,narc+2,nldr+3,nldr+4,nldr+5,narc+3,narc+4,nldr+6};

// Forme le cercle interieur
Line Loop(nlco+4) = {narc+5,narc+6};

// Crée la surface à mailler (piece de puzzle moins le cercle)
Plane Surface(nsurf+1) = {nlco+3,-(nlco+4)};

// Creation des lignes pour les conditions aux limites 
////
Physical Line(nlPh+1) = {narc+1,narc+2,nlco+4};
////
Physical Line(nlPh+2) = {nldr+1};
Physical Line(nlPh+3) = {nldr+2,nldr+3};
Physical Line(nlPh+4) = {nldr+4};
Physical Line(nlPh+5) = {nldr+5,nldr+6};

// On note 1 la surface complete
Physical Surface(nsPh+1) = {nsurf+1};

// definition de la periodicite (ouest/est)
Periodic Line {nldr+2,nldr+3}={-(nldr+6),-(nldr+5)};
// definition de la periodicite (sud/nord)
Periodic Line {nldr+1}={-(nldr+4)};

