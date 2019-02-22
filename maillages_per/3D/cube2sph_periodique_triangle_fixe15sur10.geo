pas=0.1;
// Caracteristique maillage cube 

xmin=0.;
xmax=1.;
ymin=0.;
ymax=1.;
zmin=0.;
zmax=1.;

// Caracteristique sphere aux sommets : rayon fixe
rayon_som=0.15;

/// numérotation

// Numero (moins 1) du premier point d'une arète
npar=100;
// Numero (moins 1) du premier point du cube
npcu=200;
// Numero (moins 1) du premier point de la sphère intérieure
npin=300;


// Numero (moins 1) du premier arc de cercle aux sommets
nars=500;
// Numero (moins 1) du premier arc de cercle intérieur au domaine
nari=600;

// Numero (moins 1) de la première ligne droite
nldr=700;

// Numéro (moins 1) du premier contour de face
nlf=800;
// Numéro (moins 1) du premier arc brisé aux sommets
nls=900;
// Numéro (moins 1) du premier arc brisé intérieur au domaine
nli=1000;

// Numéro (moins 1) de la première face plane
nFPl=1100;
// Numéro (moins 1) du premier huitième de sphère aux sommets : surface réglée
nSsom=1200;
// Numéro (moins 1) du premier huitième de sphère intérieur au domaine : surface réglée
nSint=1300;

// Numéro de la frontière fluide-fluide : boucle
num_surf_ff=1400;
// Numéro de la sphère intérieure : boucle
num_surf_int=1500;

// Numéro du volume à mailler
num_vol=1600;

// Numéros des ... physiques
nPhys_sol=1700;
nPhys_ff=1800;

nPhys_fluide=1900;



// Extrémités des arètes
Point(npar+1) = {xmin+rayon_som,ymin,zmin,pas};
Point(npar+2) = {xmax-rayon_som,ymin,zmin,pas};

Point(npar+3) = {xmax,ymin,zmin+rayon_som,pas};
Point(npar+4) = {xmax,ymin,zmax-rayon_som,pas};

Point(npar+5) = {xmax-rayon_som,ymin,zmax,pas};
Point(npar+6) = {xmin+rayon_som,ymin,zmax,pas};

Point(npar+7) = {xmin,ymin,zmax-rayon_som,pas};
Point(npar+8) = {xmin,ymin,zmin+rayon_som,pas};

Point(npar+9) = {xmin,ymin+rayon_som,zmin,pas};
Point(npar+10) = {xmin,ymax-rayon_som,zmin,pas};

Point(npar+11) = {xmax,ymin+rayon_som,zmin,pas};
Point(npar+12) = {xmax,ymax-rayon_som,zmin,pas};

Point(npar+13) = {xmax,ymin+rayon_som,zmax,pas};
Point(npar+14) = {xmax,ymax-rayon_som,zmax,pas};

Point(npar+15) = {xmin,ymin+rayon_som,zmax,pas};
Point(npar+16) = {xmin,ymax-rayon_som,zmax,pas};

Point(npar+17) = {xmin+rayon_som,ymax,zmin,pas};
Point(npar+18) = {xmax-rayon_som,ymax,zmin,pas};

Point(npar+19) = {xmax,ymax,zmin+rayon_som,pas};
Point(npar+20) = {xmax,ymax,zmax-rayon_som,pas};

Point(npar+21) = {xmax-rayon_som,ymax,zmax,pas};
Point(npar+22) = {xmin+rayon_som,ymax,zmax,pas};

Point(npar+23) = {xmin,ymax,zmax-rayon_som,pas};
Point(npar+24) = {xmin,ymax,zmin+rayon_som,pas};


// Sommets du cube
Point(npcu+1) = {xmin,ymin,zmin,pas};
Point(npcu+2) = {xmax,ymin,zmin,pas};
Point(npcu+3) = {xmax,ymin,zmax,pas};
Point(npcu+4) = {xmin,ymin,zmax,pas};
Point(npcu+5) = {xmin,ymax,zmin,pas};
Point(npcu+6) = {xmax,ymax,zmin,pas};
Point(npcu+7) = {xmax,ymax,zmax,pas};
Point(npcu+8) = {xmin,ymax,zmax,pas};





// Arètes : on respecte l'orientation de la base xyz dans le choix des numérotations
Line(nldr+1)={npar+1,npar+2};
Line(nldr+2)={npar+3,npar+4};
Line(nldr+3)={npar+5,npar+6};
Line(nldr+4)={npar+7,npar+8};
Line(nldr+5)={npar+9,npar+10};
Line(nldr+6)={npar+11,npar+12};
Line(nldr+7)={npar+13,npar+14};
Line(nldr+8)={npar+15,npar+16};
Line(nldr+9)={npar+17,npar+18};
Line(nldr+10)={npar+19,npar+20};
Line(nldr+11)={npar+21,npar+22};
Line(nldr+12)={npar+23,npar+24};

// Arcs centrés aux sommets : on respecte l'orientation de la base xyz dans le choix des numérotations
///ymin
Circle(nars+1)={npar+1,npcu+1,npar+8};
Circle(nars+2)={npar+3,npcu+2,npar+2};
Circle(nars+3)={npar+5,npcu+3,npar+4};
Circle(nars+4)={npar+7,npcu+4,npar+6};
///ymax
Circle(nars+5)={npar+17,npcu+5,npar+24};
Circle(nars+6)={npar+19,npcu+6,npar+18};
Circle(nars+7)={npar+21,npcu+7,npar+20};
Circle(nars+8)={npar+23,npcu+8,npar+22};
///xmin
Circle(nars+9)={npar+9,npcu+1,npar+8};
Circle(nars+10)={npar+24,npcu+5,npar+10};
Circle(nars+11)={npar+16,npcu+8,npar+23};
Circle(nars+12)={npar+7,npcu+4,npar+15};
///xmax
Circle(nars+13)={npar+11,npcu+2,npar+3};
Circle(nars+14)={npar+19,npcu+6,npar+12};
Circle(nars+15)={npar+14,npcu+7,npar+20};
Circle(nars+16)={npar+4,npcu+3,npar+13};
///zmin
Circle(nars+17)={npar+1,npcu+1,npar+9};
Circle(nars+18)={npar+11,npcu+2,npar+2};
Circle(nars+19)={npar+18,npcu+6,npar+12};
Circle(nars+20)={npar+10,npcu+5,npar+17};
///zmax
Circle(nars+21)={npar+6,npcu+4,npar+15};
Circle(nars+22)={npar+13,npcu+3,npar+5};
Circle(nars+23)={npar+21,npcu+7,npar+14};
Circle(nars+24)={npar+16,npcu+8,npar+22};






// Contours des faces puis faces

//xmin : 1
Line Loop(nlf+1)={nldr+5,-(nars+10),-(nldr+12),-(nars+11),-(nldr+8),-(nars+12),nldr+4,-(nars+9)};
Plane Surface(nFPl+1)={nlf+1};

//xmax : 6
Line Loop(nlf+6)={nldr+6,-(nars+14),nldr+10,-(nars+15),-(nldr+7),-(nars+16),-(nldr+2),-(nars+13)};
Plane Surface(nFPl+6)={nlf+6};

//ymin : 3
Line Loop(nlf+3)={nldr+1,-(nars+2),nldr+2,-(nars+3),nldr+3,-(nars+4),nldr+4,-(nars+1)};
Plane Surface(nFPl+3)={nlf+3};

//ymax : 4
Line Loop(nlf+4)={nldr+9,-(nars+6),nldr+10,-(nars+7),nldr+11,-(nars+8),nldr+12,-(nars+5)};
Plane Surface(nFPl+4)={nlf+4};

//zmin : 2
Line Loop(nlf+2)={nldr+1,-(nars+18),nldr+6,-(nars+19),-(nldr+9),-(nars+20),-(nldr+5),-(nars+17)};
Plane Surface(nFPl+2)={nlf+2};

//zmax : 5
Line Loop(nlf+5)={-(nldr+3),-(nars+22),nldr+7,-(nars+23),nldr+11,-(nars+24),-(nldr+8),-(nars+21)};
Plane Surface(nFPl+5)={nlf+5};



// Huitièmes des sphères centrés aux sommets : contours puis surfaces
//y=ymin
Line Loop(nls+1)={nars+17,nars+9,-(nars+1)};
Line Loop(nls+2)={nars+18,-(nars+2),-(nars+13)};
Line Loop(nls+3)={nars+22,nars+3,nars+16};
Line Loop(nls+4)={nars+21,-(nars+12),nars+4};
//y=ymax
Line Loop(nls+5)={nars+20,nars+5,nars+10};
Line Loop(nls+6)={nars+19,-(nars+14),nars+6};
Line Loop(nls+7)={(nars+23),nars+15,-(nars+7)};
Line Loop(nls+8)={(nars+24),-(nars+8),-(nars+11)};

Ruled Surface(nSsom+1)={nls+1};
Ruled Surface(nSsom+2)={nls+2};
Ruled Surface(nSsom+3)={nls+3};
Ruled Surface(nSsom+4)={nls+4};
Ruled Surface(nSsom+5)={nls+5};
Ruled Surface(nSsom+6)={nls+6};
Ruled Surface(nSsom+7)={nls+7};
Ruled Surface(nSsom+8)={nls+8};



// Surface fluide-fluide : assemblage des faces planes et périodicité

////
Surface Loop(num_surf_ff)={(nFPl+1),(nFPl+5),-(nFPl+3),(nFPl+6),(nFPl+2),-(nFPl+4),nSsom+1,nSsom+2,nSsom+3,nSsom+4,nSsom+5,nSsom+6,nSsom+7,nSsom+8};
////

/// ouest-est : 1-6
Periodic Surface nFPl+1 {(nldr+5),(nars+12),-(nldr+12),(nars+9),(nldr+8),(nars+11),(nldr+4),(nars+10)}= nFPl+6 {(nldr+6),(nars+16),(nldr+10),(nars+13),(nldr+7),(nars+15),-(nldr+2),(nars+14)};

/// sud-nord : 3-4
Periodic Surface nFPl+3 {(nldr+1),(nars+1),(nldr+2),(nars+2),(nldr+3),(nars+3),(nldr+4),(nars+4)}= nFPl+4 {(nldr+9),(nars+5),(nldr+10),(nars+6),(nldr+11),(nars+7),(nldr+12),(nars+8)};

/// bas-haut : 2-5
Periodic Surface nFPl+2 {(nldr+1),(nars+17),(nldr+6),(nars+20),(nldr+9),(nars+19),(nldr+5),(nars+18)}= nFPl+5 {-(nldr+3),(nars+21),(nldr+7),(nars+24),-(nldr+11),(nars+23),(nldr+8),(nars+22)};



// Volume à mailler

Volume(num_vol)={num_surf_ff};//nSsom+1,nSsom+2,nSsom+3,nSsom+4,nSsom+5,nSsom+6,nSsom+7,nSsom+8,(nFPl+1),(nFPl+3),(nFPl+5),(nFPl+6),(nFPl+4),(nFPl+2)};//

// Surfaces et volume physiques

Physical Surface(nPhys_sol)={nSsom+1,nSsom+2,nSsom+3,nSsom+4,nSsom+5,nSsom+6,nSsom+7,nSsom+8};
Physical Surface(nPhys_ff)={(nFPl+1),(nFPl+3),(nFPl+5),(nFPl+6),(nFPl+4),(nFPl+2)};//,num_surf_ff

Physical Volume(nPhys_fluide)={num_vol};






















