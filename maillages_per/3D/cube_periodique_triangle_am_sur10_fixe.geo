// Caracteristique maillage cube 
pas_cube=0.1;// pas_cube : pas pour mailler le cube
xmin=0.;
xmax=1.;
ymin=0.;
ymax=1.;
zmin=0.;
zmax=1.;

// Numero (moins 1) du premier point définissant le cube
npc=100;
// Numero (moins 1) de la première ligne définissant le cube
nlc=200;
// Numero (moins 1) de la première Surface plane définissant le cube
nPSc=300;
// Numero (moins 1) de la première Surface définissant le cube
nsc=400;
// Numero de la surface totale définissant le cube
num_surf_loop_cube=600;
// Numero du volume
num_vol=500;
// Numero des surfaces physiques
num_phys=10;

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
// Début du maillage du cube 
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
// Creation des 8 points du bord du cube
Point(npc+1) = {xmin,ymin,zmin,pas_cube};
Point(npc+2) = {xmax,ymin,zmin,pas_cube};
Point(npc+3) = {xmax,ymax,zmin,pas_cube};
Point(npc+4) = {xmin,ymax,zmin,pas_cube};
Point(npc+5) = {xmin,ymin,zmax,pas_cube};
Point(npc+6) = {xmax,ymin,zmax,pas_cube};
Point(npc+7) = {xmax,ymax,zmax,pas_cube};
Point(npc+8) = {xmin,ymax,zmax,pas_cube};

// Creation des 8 arêtes du cube
Line(nlc+1) = {npc+4,npc+3};
Line(nlc+2) = {npc+3,npc+2};
Line(nlc+3) = {npc+2,npc+1};
Line(nlc+4) = {npc+1,npc+4};
Line(nlc+6) = {npc+5,npc+6};
Line(nlc+7) = {npc+6,npc+7};
Line(nlc+8) = {npc+7,npc+8};
Line(nlc+9) = {npc+8,npc+5};
Line(nlc+10) = {npc+1,npc+5};
Line(nlc+11) = {npc+4,npc+8};
Line(nlc+12) = {npc+2,npc+6};
Line(nlc+13) = {npc+3,npc+7};

// Definition du contour fermé et de la surface située en zmin
Line Loop(nlc+14) = {nlc+3,nlc+4,nlc+1,nlc+2};
Plane Surface(nPSc+15) = {nlc+14};
// Definition du contour fermé et de la surface située en zmax
Line Loop(nlc+16) = {nlc+6,nlc+7,nlc+8,nlc+9};
Plane Surface(nPSc+17) = {nlc+16};
// Definition du contour fermé et de la surface située en xmin
Line Loop(nlc+18) = {nlc+10,-(nlc+9),-(nlc+11),-(nlc+4)};
Plane Surface(nPSc+19) = {nlc+18};
// Definition du contour fermé et de la surface située en ymax
Line Loop(nlc+20) = {nlc+8,-(nlc+11),nlc+1,nlc+13};
Plane Surface(nPSc+21) = {nlc+20};
// Definition du contour fermé et de la surface située en xmax
Line Loop(nlc+22) = {nlc+12,nlc+7,-(nlc+13),nlc+2};
Plane Surface(nPSc+23) = {nlc+22};
// Definition du contour fermé et de la surface située en ymin
Line Loop(nlc+24) = {nlc+6,-(nlc+12),nlc+3,nlc+10};
Plane Surface(nPSc+25) = {nlc+24};

// // On impose la periodicité entre les surfaces d'équations xmin et xmax
Periodic Surface nPSc+23 {nlc+12,-(nlc+7), nlc+13, -(nlc+2)} = nPSc+19 {nlc+10, nlc+9, nlc+11, nlc+4};
// On impose la periodicité entre les surfaces d'équations ymin et ymax
Periodic Surface nPSc+21 {nlc+11,-(nlc+8), nlc+13, -(nlc+1)} = nPSc+25 {nlc+10, nlc+6, nlc+12, nlc+3};
// // On impose la periodicité entre les surfaces d'équations zmin et zmax
Periodic Surface nPSc+17 {-(nlc+8),-(nlc+7), -(nlc+6), -(nlc+9)} = nPSc+15 {nlc+1, nlc+2, nlc+3, nlc+4};


Surface Loop(num_surf_loop_cube) = {nPSc+17,-(nPSc+25),-(nPSc+23),-(nPSc+21),nPSc+19,nPSc+15};

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
// Fin du maillage du cube 
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////





/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
// Creation du volume à mailler
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
Volume(num_vol) = {num_surf_loop_cube};


/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
// Definition des surfaces physiques
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
Physical Surface(num_phys+1) = {nPSc+15}; //surface située en zmin
Physical Surface(num_phys+2) = {nPSc+17}; //surface située en zmax
Physical Surface(num_phys+3) = {nPSc+19}; // surface située en xmin
Physical Surface(num_phys+4) = {nPSc+21}; // surface située en ymax
Physical Surface(num_phys+5) = {nPSc+23}; // surface située en xmax
Physical Surface(num_phys+6) = {nPSc+25}; // surface située en ymin

Physical Volume(num_phys+7) = {num_vol}; 
