### Une fonction, pour mailler un domaine avec une inclusion circulaire. Il en existe d'autres, voir homog_pod_multi ###

def raffinemment_maillage(r,mesh):
    markers = CellFunction("bool", mesh)
    markers.set_all(False)
    for c in cells(mesh):
	# Mark cells with facet midpoints near y == 1.0
	for f in facets(c):
	    if (sqrt((f.midpoint()[0]-c_x)**2+(f.midpoint()[1]-c_y)**2)<=1.2*r):
		markers[c] = True
    #new_mesh=refine(mesh, markers, redistribute=False)  
    new_mesh=refine(mesh, markers, redistribute=True)  
    
    return new_mesh

### Classes pour caractériser les bordures Gamma_sf et Gamma_ff , cas circulaire ###

# Sub domain for Periodic boundary condition
class PeriodicBoundary(SubDomain):
        
        # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        return on_boundary and (near(x[0],xinf,tol) or near(x[1],yinf,tol))

            
            # Map right boundary (H) to left boundary (G)
    def map(self, x, y):
        if (near(x[0],xsup,tol)):
            y[0] = x[0] - 1.0
            y[1] = x[1]
                
        else :
            y[0]=x[0]
            y[1] = x[1] - 1.0

# A renommer : interface #

class Obstacle(SubDomain):
 def inside(self, x, on_boundary):
  return (on_boundary and between((x[0]-c_x), (-r-tol, r+tol)) and between((x[1]-c_x), (-r-tol, r+tol)))

## Cas géométrique plus général, à réécrire ##

# Domaine d'intégration de chi, éventuellement pour des géométries plus compliquées que circulaire #

class Omega_fluide(SubDomain):
 def inside(self, x, on_boundary):
  return (((sqrt((x[0] - c_x)**2 + (x[1] - c_y)**2)-r) >= 0.0))


