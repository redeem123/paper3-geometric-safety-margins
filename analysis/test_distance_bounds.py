import numpy as np
from scipy.optimize import minimize_scalar, differential_evolution
from distance_bounds import Terrain,distance,box_bounds,patch_value,split_patch

def terrain(H,x=[0.,10.],y=[0.,10.]):
 X,Y=np.meshgrid(x,y);return Terrain(dict(terrainX=X,terrainY=Y,terrainGrid=np.asarray(H,float)))
# Independent scalar minimization of convex point-to-box distance along a segment.
rng=np.random.default_rng(20261112)
for _ in range(500):
 P=rng.uniform(-20,20,(2,3));lo=rng.uniform(-10,10,3);hi=lo+rng.uniform(.01,10,3)
 lb,q=box_bounds(P,lo,hi)
 def f(t):
  p=P[0]+t*(P[1]-P[0]);return np.sum((p-np.clip(p,lo,hi))**2)
 rr=minimize_scalar(f,bounds=(0,1),method='bounded',options={'xatol':1e-14})
 expected=np.sqrt(min(f(0),f(1),rr.fun));assert abs((lb[0]+1e-9)-expected)<1e-6
flat=terrain([[0,0],[0,0]])
r=distance(np.array([[1,1,8],[9,9,8.]]),flat,tolerance=1e-5)
assert r['lower']<=8<=r['upper'] and r['gap']<=1e-5
# Plane z=x: the vertical gap10 gives perpendicular distance10/sqrt2.
slope=terrain([[0,10],[0,10]])
r=distance(np.array([[1,1,11],[5,9,15.]]),slope,tolerance=1e-4)
expected=10/np.sqrt(2);assert r['lower']<=expected<=r['upper'] and r['converged']
# Vertical route over a bilinear saddle, compared with an independent global search.
T=terrain([[0,6],[7,1]])
P=np.array([[2.,3.,12.],[7.,8.,14.]])
r=distance(P,T,tolerance=1e-4)
p=T.patch(0)
def f(w):
 u,v,t=w;point=np.array([10*u,10*v,patch_value(p,[[10*u,10*v]])[0]])
 return np.linalg.norm(point-(P[0]+t*(P[1]-P[0])))
opt=differential_evolution(f,[(0,1)]*3,tol=1e-11,seed=230,polish=True)
assert r['lower']<=opt.fun+1e-8 and r['upper']>=opt.fun-1e-8 and r['converged']
# Subdivision reproduces the same bilinear surface.
for child in split_patch(p):
 x0,x1,y0,y1,*_=child;xy=np.column_stack([rng.uniform(x0,x1,50),rng.uniform(y0,y1,50)])
 assert np.max(abs(patch_value(child,xy)-patch_value(p,xy)))<1e-12
print('PASS:500 segment/AABB cases, flat and sloped planes, independent bilinear global minimization, exact patch subdivision.')
# Parallel paths whose closest points are not their endpoints, and intersection.
r=distance(np.array([[1.,5.,5.],[9.,5.,5.]]),flat,tolerance=1e-4)
assert r['lower']<=5<=r['upper'] and r['converged']
r=distance(np.array([[5.,5.,-2.],[5.,5.,2.]]),flat,tolerance=1e-4)
assert r['lower']==0 and r['upper']<1e-7
from triangle_bound import patch_triangle_bound,seg_edges,polyline_triangle
dd,_=polyline_triangle(np.array([[-2.,5.,5.],[12.,5.,5.]]),np.array([0.,0.,0.]),np.array([10.,0.,0.]),np.array([0.,10.,0.]))
assert abs(dd-5)<1e-12
p=np.array([[0.,-5e-9,0.],[10.,5e-9,0.]])
ds,qq=seg_edges(p,np.array([0.,0.,0.]),np.array([10.,0.,0.]))
assert ds[0]<1e-20
for _ in range(5):
 T=terrain(rng.uniform(-3,8,(2,2)));p=T.patch(0);P=rng.uniform([0,0,10],[10,10,18],(2,3))
 r=distance(P,T,tolerance=1e-4)
 def fun(w):
  u,v,t=w;z=patch_value(p,[[10*u,10*v]])[0]
  return np.linalg.norm(np.array([10*u,10*v,z])-(P[0]+t*(P[1]-P[0])))
 o=differential_evolution(fun,[(0,1)]*3,tol=1e-10,seed=75,polish=True)
 assert r['lower']<=o.fun+1e-7 and r['upper']>=o.fun-1e-7 and r['converged'],(r,o.fun)
print('PASS: interior parallel minima, intersections, near-parallel edges and5 additional independent bilinear minimizations.')
