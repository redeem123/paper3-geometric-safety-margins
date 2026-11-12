from pathlib import Path
import sys,json,heapq,csv,time
import numpy as np
from scipy.io import loadmat
from scipy.interpolate import RegularGridInterpolator
W=Path(__file__).resolve().parent;M=W;O=W/'recomputed_geometry';O.mkdir(exist_ok=True)
from portable_routes import read_original_route
from distance_bounds import Terrain,box_bounds
from triangle_bound import polyline_triangle
from native_geometry import init,distance
rows=list(csv.DictReader((M/'verified_runs.csv').open()));original={(int(r['scenario']),r['method'],int(r['run'])):r for r in csv.DictReader((M/'euclidean/route_distances.csv').open())}
def trimin(P,T,opposite=False):
 heap=[(0.,0)];upper=float('inf');cnt=0
 while heap:
  lb,node=heapq.heappop(heap)
  if lb>=upper:continue
  lo,hi,children,_=T.nodes[node]
  if children:
   b,_=box_bounds(P,np.array([T.nodes[n][0] for n in children]),np.array([T.nodes[n][1] for n in children]))
   for n,d in zip(children,b):
    if d<upper:heapq.heappush(heap,(float(d),n))
  else:
   x0,x1,y0,y1,z00,z10,z01,z11=T.patch(node);a=np.array([x0,y0,z00]);b=np.array([x1,y0,z10]);c=np.array([x1,y1,z11]);d=np.array([x0,y1,z01])
   faces=[(a,b,d),(b,c,d)] if opposite else [(a,b,c),(a,c,d)]
   upper=min(upper,*(polyline_triangle(P,*f)[0] for f in faces));cnt+=1
 return upper,cnt
# Analytic controls and near-threshold/edge cases for the actual native engine.
tests=[]
for slope in [0,.2,1,10]:
 x=np.array([0.,100.]);y=x;H=np.array([[0,100*slope],[0,100*slope]]);init(x,y,H)
 for delta in [-1e-5,0,1e-5]:
  h=(8+delta)*np.sqrt(1+slope*slope);P=np.array([[1,20,slope+h],[2,80,2*slope+h]]);r=distance(P,1e-7);exact=8+delta
  assert r[0]<=exact+1e-8 and r[1]>=exact-1e-8;tests.append(dict(slope=slope,exact=exact,lower=r[0],upper=r[1]))
# Original single cells subdivided exactly, including boundaries and saddle.
rng=np.random.default_rng(1801)
for k in range(10):
 x=np.array([0.,5.,10.]);X,Y=np.meshgrid(x,x);H=(X*Y*.03+X*.2-Y*.1);P=np.array([[0.,0.,9.],[10.,10.,15.]])+[0,0,k*.01];init(x,x,H);a=distance(P,1e-5)
 xx=np.linspace(0,10,9);XX,YY=np.meshgrid(xx,xx);HH=XX*YY*.03+XX*.2-YY*.1;init(xx,xx,HH);b=distance(P,1e-5)
 assert a[0]<=b[1] and b[0]<=a[1]
(O/'analytic_geometry_tests.json').write_text(json.dumps(dict(passed=True,nearThresholdPlanes=tests,cellBoundarySubdivisionChecks=10),indent=2))
alt=[];profiles=[]
for s in [101,102,103]:
 group=[r for r in rows if int(r['scenario'])==s];env=read_original_route(s,group[0]['method'],int(group[0]['run']))['env'];T=Terrain(env);x=T.x;y=T.y;H=T.H
 # Same vertices, two triangulation orientations. No terrain heights or paths changed.
 for rr in group:
  P=read_original_route(s,rr['method'],int(rr['run']))['path'];v,c=trimin(P,T);v2,c2=trimin(P,T,True);o=original[s,rr['method'],int(rr['run'])]
  alt.append(dict(scenario=s,method=rr['method'],run=int(rr['run']),bilinearLower=float(o['lower']),bilinearUpper=float(o['upper']),triangleMain=v,triangleOther=v2))
 print('triangles',s,flush=True)
 # Regridding keeps the original surface for refinement, and coarsens the same target for sensitivity.
 selected=[r for r in group if int(r['run'])==1]
 for factor in [2,.5]:
  if factor==2:xx=np.sort(np.r_[x,(x[1:]+x[:-1])/2]);yy=np.sort(np.r_[y,(y[1:]+y[:-1])/2])
  else:xx=x[np.r_[np.arange(0,len(x)-1,2),len(x)-1]];yy=y[np.r_[np.arange(0,len(y)-1,2),len(y)-1]]
  XX,YY=np.meshgrid(xx,yy);HH=T.interp(np.stack([YY,XX],axis=-1));ip=RegularGridInterpolator((yy,xx),HH)
  X0,Y0=np.meshgrid(x,y);heightError=float(np.max(abs(ip(np.stack([Y0,X0],axis=-1))-H)))
  init(xx,yy,HH)
  for rr in selected:
   P=read_original_route(s,rr['method'],int(rr['run']))['path'];d=distance(P,1e-4);o=original[s,rr['method'],1]
   if factor==2:assert d[0]<=float(o['upper'])+1e-7 and float(o['lower'])<=d[1]+1e-7
   profiles.append(dict(kind='resolution',scenario=s,method=rr['method'],run=1,factor=factor,grid=len(xx),heightError=heightError,lower=d[0],upper=d[1]))
 # Spatial exposure as bounded arc length of {d(s)<8};distance to surface is1-Lipschitz.
 init(x,y,H)
 for rr in selected:
  P=read_original_route(s,rr['method'],int(rr['run']))['path'];L=float(np.linalg.norm(np.diff(P,axis=0),axis=1).sum());heap=[];lowerLength=0.;unknown=0.;serial=0;calls=0
  def add(a,b):
   global lowerLength,unknown,serial,calls
   length=float(np.linalg.norm(b-a));d=distance([(a+b)/2],1e-4);calls+=1;radius=length/2
   if d[1]+radius<8:lowerLength+=length
   elif d[0]-radius>=8:pass
   else:unknown+=length;serial+=1;heapq.heappush(heap,(-length,serial,a,b))
  for a,b in zip(P[:-1],P[1:]):add(a,b)
  while unknown/L>1e-3:
   length,_,a,b=heapq.heappop(heap);unknown+=length;mid=(a+b)/2;add(a,mid);add(mid,b)
   assert calls<100000
  profiles.append(dict(kind='spatial_profile',scenario=s,method=rr['method'],run=1,lower=max(0,lowerLength/L-1e-12),upper=min(1,max(lowerLength,lowerLength+unknown)/L+1e-12),pointQueries=calls))
 print('profiles',s,flush=True)
with (O/'surface_alternatives.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(alt[0]));w.writeheader();w.writerows(alt)
(O/'spatial_profiles_resolution.json').write_text(json.dumps(profiles,indent=2))
print('Completed630 x2 triangulations,42 regridding queries,21 spatial profiles.',flush=True)
