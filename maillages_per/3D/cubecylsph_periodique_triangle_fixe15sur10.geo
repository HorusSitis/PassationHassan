pas=0.1;
// Caracteristique maillage cube 

xmin=0.;
xmax=1.;
ymin=0.;
ymax=1.;
zmin=0.;
zmax=1.;

// Caracteristique sphere interieure
xc=0.5;
yc=0.5;
zc=0.5;

rayon_sph=0.15;

// Caracteristique sphere aux sommets : rayon fixe
rayon_cyl=0;

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
// Numéro (moins 1) du premier contour de quart de cylindre
nlc=900;
// Numéro (moins 1) du premier arc brisé intérieur au domaine
nli=1000;

// Numéro (moins 1) de la première face plane
nFPl=1100;
// Numéro (moins 1) du premier quart de cylindre : surface réglée
nSsom=1200;
// Numéro (moins 1) du premier huitième de sphère intérieur au domaine : surface réglée
nSint=1300;

// Numéro de la frontière fluide-fluide : boucle
num_surf_ext=1400;
// Numéro de la sphère intérieure : boucle
num_surf_int=1500;

// Numéro du volume à mailler
num_vol=1600;

// Numéros des ... physiques
nPhys_sol=1700;
nPhys_ff=1800;

nPhys_fluide=1900;



// Extrémités des arètes
Point(npar+1)={xmin,ymin,zmin,pas};
Point(npar+2)={xmin,ymin,zmin,pas};
Point(npar+3)={xmax,ymin,zmin,pas};
Point(npar+4)={xmax,ymin,zmin,pas};
Point(npar+5)={xmax,ymin,zmax,pas};
Point(npar+6)={xmax,ymin,zmax,pas};
Point(npar+7)={xmin,ymin,zmax,pas};
Point(npar+8)={xmin,ymin,zmax,pas};
Point(npar+9)={xmin,ymax,zmin,pas};
Point(npar+10)={xmin,ymax,zmin,pas};
Point(npar+11)={xmax,ymax,zmin,pas};
Point(npar+12)={xmax,ymax,zmin,pas};
Point(npar+13)={xmax,ymax,zmax,pas};
Point(npar+14)={xmax,ymax,zmax,pas};
Point(npar+15)={xmin,ymax,zmax,pas};
Point(npar+16)={xmin,ymax,zmax,pas};



// Sommets du cube
Point(npcu+1) = {xmin,ymin,zmin,pas};
Point(npcu+2) = {xmax,ymin,zmin,pas};
Point(npcu+3) = {xmax,ymin,zmax,pas};
Point(npcu+4) = {xmin,ymin,zmax,pas};
Point(npcu+5) = {xmin,ymax,zmin,pas};
Point(npcu+6) = {xmax,ymax,zmin,pas};
Point(npcu+7) = {xmax,ymax,zmax,pas};
Point(npcu+8) = {xmin,ymax,zmax,pas};


// Points définissant la sphère intérieure
Point(npin+1) = {xc,yc,zc,pas};
Point(npin+2) = {xc+rayon_sph,yc,zc,pas};
Point(npin+3) = {xc,yc+rayon_sph,zc,pas};
Point(npin+4) = {xc-rayon_sph,yc,zc,pas};
Point(npin+5) = {xc,yc-rayon_sph,zc,pas};
Point(npin+6) = {xc,yc,zc+rayon_sph,pas};
Point(npin+7) = {xc,yc,zc-rayon_sph,pas};


// Arètes : on respecte l'orientation de la base xyz dans le choix des numérotations
Line(nldr+1)={npar+1,npar+3};
Line(nldr+2)={npar+3,npar+5};
Line(nldr+3)={npar+5,npar+7};
Line(nldr+4)={npar+7,npar+1};

Line(nldr+5)={npar+1,npar+9};
Line(nldr+6)={npar+1,npar+9};
Line(nldr+7)={npar+3,npar+11};
Line(nldr+8)={npar+3,npar+11};
Line(nldr+9)={npar+5,npar+13};
Line(nldr+10)={npar+5,npar+13};
Line(nldr+11)={npar+7,npar+15};
Line(nldr+12)={npar+7,npar+15};

Line(nldr+13)={npar+9,npar+11};
Line(nldr+14)={npar+11,npar+13};
Line(nldr+15)={npar+13,npar+15};
Line(nldr+16)={npar+15,npar+9};







// Arcs centrés à l'intérieur du domaine

///z=
Circle(nari+1)={npin+2,npin+1,npin+3};
Circle(nari+2)={npin+3,npin+1,npin+4};
Circle(nari+3)={npin+4,npin+1,npin+5};
Circle(nari+4)={npin+5,npin+1,npin+2};
///y=
Circle(nari+5)={npin+2,npin+1,npin+6};
Circle(nari+6)={npin+6,npin+1,npin+4};
Circle(nari+7)={npin+4,npin+1,npin+7};
Circle(nari+8)={npin+7,npin+1,npin+2};
///x=
Circle(nari+9)={npin+6,npin+1,npin+5};
Circle(nari+10)={npin+5,npin+1,npin+7};
Circle(nari+11)={npin+7,npin+1,npin+3};
Circle(nari+12)={npin+3,npin+1,npin+6};


// Contours des faces puis faces : orientation vers l'extérieur pour xzminmax et ymax, inversement pour yminmax

//xmin : 1
Line Loop(nlf+1)={nldr+16,-(nldr+5),-(nldr+4),nldr+12};
Plane Surface(nFPl+1)={nlf+1};

//xmax : 6
Line Loop(nlf+6)={nldr+14,-(nldr+9),-(nldr+2),nldr+8};
Plane Surface(nFPl+6)={nlf+6};

//ymin : 3
Line Loop(nlf+3)={-(nldr+1),-(nldr+4),-(nldr+3),-(nldr+2)};
Plane Surface(nFPl+3)={nlf+3};

//ymax : 4
Line Loop(nlf+4)={nldr+13,nldr+14,nldr+15,nldr+16};
Plane Surface(nFPl+4)={nlf+4};

//zmin : 2
Line Loop(nlf+2)={nldr+1,nldr+7,-(nldr+13),-(nldr+6)};
Plane Surface(nFPl+2)={nlf+2};

//zmax : 5
Line Loop(nlf+5)={-(nldr+3),nldr+10,nldr+15,-(nldr+11)};
Plane Surface(nFPl+5)={nlf+5};






// Sphère intérieure
//z>=zc
Line Loop(nli+1)={nari+1,nari+12,-(nari+5)};
Line Loop(nli+2)={nari+2,-(nari+6),-(nari+12)};
Line Loop(nli+3)={nari+3,-(nari+9),nari+6};
Line Loop(nli+4)={nari+4,nari+5,nari+9};
//z<=zc
Line Loop(nli+5)={nari+4,-(nari+8),-(nari+10)};
Line Loop(nli+6)={nari+3,nari+10,-(nari+7)};
Line Loop(nli+7)={nari+2,nari+7,nari+11};
Line Loop(nli+8)={nari+1,-(nari+11),nari+8};

Ruled Surface(nSint+1)={nli+1};
Ruled Surface(nSint+2)={nli+2};
Ruled Surface(nSint+3)={nli+3};
Ruled Surface(nSint+4)={nli+4};
Ruled Surface(nSint+5)={nli+5};
Ruled Surface(nSint+6)={nli+6};
Ruled Surface(nSint+7)={nli+7};
Ruled Surface(nSint+8)={nli+8};

Surface Loop(num_surf_int)={nSint+1,-(nSint+2),nSint+3,-(nSint+4),-(nSint+5),(nSint+6),-(nSint+7),(nSint+8)};

// Frontière extérieure du domaine fluide : assemblage des faces planes et quarts de cylindre, périodicité

////
Surface Loop(num_surf_ext)={(nFPl+1),(nFPl+2),-(nFPl+3),-(nFPl+6),-(nFPl+5),(nFPl+4)};
////

/// ouest-est : 1-6
Periodic Surface nFPl+1 {(nldr+5),(nldr+16),(nldr+12),(nldr+4)}= nFPl+6 {(nldr+8),-(nldr+14),(nldr+9),-(nldr+2)};

/// sud-nord : 3-4
Periodic Surface nFPl+3 {(nldr+1),(nldr+2),(nldr+3),(nldr+4)}=nFPl+4 {(nldr+13),(nldr+14),(nldr+15),(nldr+16)};

/// bas-haut : 2-5
Periodic Surface nFPl+2 {(nldr+1),(nldr+7),(nldr+13),(nldr+6)}= nFPl+5 {-(nldr+3),(nldr+10),-(nldr+15),(nldr+11)};



// Volume à mailler

Volume(num_vol)={num_surf_int,num_surf_ext};

// Surfaces et volume physiques

Physical Surface(nPhys_sol)={num_surf_int};
Physical Surface(nPhys_ff)={(nFPl+1),(nFPl+3),(nFPl+5),(nFPl+6),(nFPl+4),(nFPl+2)};

Physical Volume(nPhys_fluide)={num_vol};






















