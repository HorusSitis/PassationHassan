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

rayon_int=0.05;

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


// Points définissant la sphère intérieure
Point(npin+1) = {xc,yc,zc,pas};
Point(npin+2) = {xc+rayon_int,yc,zc,pas};
Point(npin+3) = {xc,yc+rayon_int,zc,pas};
Point(npin+4) = {xc-rayon_int,yc,zc,pas};
Point(npin+5) = {xc,yc-rayon_int,zc,pas};
Point(npin+6) = {xc,yc,zc+rayon_int,pas};
Point(npin+7) = {xc,yc,zc-rayon_int,pas};


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
Circle(nars+23)={npar+14,npcu+7,npar+22};
Circle(nars+24)={npar+16,npcu+8,npar+21};



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





// Contours des faces puis faces






// Huitièmes des sphères centrés aux sommets : contours puis surfaces




// Sphère intérieure





























