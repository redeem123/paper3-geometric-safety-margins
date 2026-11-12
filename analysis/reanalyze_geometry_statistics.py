from pathlib import Path
import sys,json,csv,hashlib,itertools,heapq,time
import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.stats import mannwhitneyu,spearmanr
W=Path(__file__).resolve().parent;M=W;O=W/'recomputed_analysis';O.mkdir(exist_ok=True)
from portable_routes import read_original_route
from clearance_profiles import from_route
from distance_bounds import Terrain,distance,box_bounds
A=pd.read_csv(M/'verified_runs.csv');D=pd.read_csv(M/'euclidean/route_distances.csv');A=A.merge(D[['scenario','method','run','lower','upper','seconds']],on=['scenario','method','run'],suffixes=('','_assessment'))
methods=list(A.method.drop_duplicates());maps=[101,102,103]
def wilson(k,n):
 z=1.959963984540054;p=k/n;c=(p+z*z/(2*n))/(1+z*z/n);half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n);return [max(0,c-half),min(1,c+half)]
def clip(P,alpha):
 ell=np.linalg.norm(np.diff(P,axis=0),axis=1);s=np.r_[0,np.cumsum(ell)];a=alpha*s[-1];b=(1-alpha)*s[-1]
 q=np.array([np.interp(a,s,P[:,j]) for j in range(3)]);r=np.array([np.interp(b,s,P[:,j]) for j in range(3)])
 return np.vstack([q,P[(s>a)&(s<b)],r])
rows=[];en=[];terrain={};per=[];curves=[]
for row in A.itertuples():
 a=read_original_route(row.scenario,row.method,row.run);P=a['path'];e=a['env'];s=row.scenario
 if s not in terrain:
  terrain[s]=Terrain(e);T=terrain[s]
  for label,p in [('start',P[0]),('goal',P[-1])]:
   # Equivalent zero-length-point query via a tiny upwards segment brackets
   # point distance within its length; account explicitly for that length.
   pp=np.vstack([p,p+[0,0,1e-7]]);b=distance(pp,T,tolerance=1e-5)
   en.append(dict(scenario=s,endpoint=label,lower=b['lower'],upper=b['upper']+1e-7))
 T=terrain[s];j=json.loads((M/f'euclidean/results/scenario{s}_{row.method}_run{row.run}.json').read_text())
 q=np.array(j['surfacePoint']);rp=np.array(j['pathPoint']);ix=np.clip(np.searchsorted(T.x,q[0],side='right')-1,0,len(T.x)-2);iy=np.clip(np.searchsorted(T.y,q[1],side='right')-1,0,len(T.y)-2)
 u=(q[0]-T.x[ix])/(T.x[ix+1]-T.x[ix]);v=(q[1]-T.y[iy])/(T.y[iy+1]-T.y[iy]);h=T.H[iy:iy+2,ix:ix+2]
 gx=((h[0,1]-h[0,0])*(1-v)+(h[1,1]-h[1,0])*v)/(T.x[ix+1]-T.x[ix]);gy=((h[1,0]-h[0,0])*(1-u)+(h[1,1]-h[0,1])*u)/(T.y[iy+1]-T.y[iy]);hw=rp[2]-T.z(rp[:2])[0]
 out=dict(scenario=s,method=row.method,run=row.run,verticalMin=row.clearance,lower=row.lower,upper=row.upper,ratioLow=row.clearance/row.upper,ratioHigh=row.clearance/row.lower,witnessVertical=hw,slopeNorm=np.hypot(gx,gy))
 for alpha in [.01,.05,.1]:out[f'trim{int(100*alpha)}']=from_route(clip(P,alpha),e).minimum
 rows.append(out)
 if row.run==1:print('analyzed',s,row.method,flush=True)
 profile=from_route(P,e);cc=np.linspace(8,80,145);curve=profile.cdf(cc)
 curves.extend([dict(scenario=s,method=row.method,run=row.run,c=c,F=f) for c,f in zip(cc,curve)])
pd.DataFrame(rows).to_csv(O/'geometry_profiles.csv',index=False);pd.DataFrame(en).to_csv(O/'endpoint_distances.csv',index=False)
# All unordered tests use independent run samples, not arbitrary matching by run index.
for s in maps:
 for a,b in itertools.combinations(methods,2):
  x=A[(A.scenario==s)&(A.method==a)].J.to_numpy();y=A[(A.scenario==s)&(A.method==b)].J.to_numpy();u,p=mannwhitneyu(x,y,alternative='two-sided')
  per.append(dict(scenario=s,A=a,B=b,medianA=np.median(x),medianB=np.median(y),medianDifference=np.median(x)-np.median(y),probabilityALower=np.mean(x[:,None]<y)+.5*np.mean(x[:,None]==y),cliffsDelta=2*u/900-1,p=p))
order=np.argsort([r['p'] for r in per]);prev=0
for rank,i in enumerate(order):prev=max(prev,min(1,per[i]['p']*(len(per)-rank)));per[i]['holmP']=prev
pd.DataFrame(per).to_csv(O/'all_pair_tests.csv',index=False)
# Relative materiality means both fR and eta differ by at least delta, using symmetric scales.
rev=[];examples=[]
for s in maps:
 for a,b in itertools.combinations(methods,2):
  x=A[(A.scenario==s)&(A.method==a)].copy();y=A[(A.scenario==s)&(A.method==b)].copy();fx=x.fR.to_numpy()[:,None];fy=y.fR.to_numpy()[None,:];ex=(x.fR/x.fL).to_numpy()[:,None];ey=(y.fR/y.fL).to_numpy()[None,:];diff=(fx-fy)*(ex-ey)<0;ties=(abs(fx-fy)<=1e-7)|(abs(ex-ey)<=1e-7)
  out=dict(scenario=s,A=a,B=b)
  for delta in [0,.01,.05,.1]:
   good=diff&~ties&(2*abs(fx-fy)/(fx+fy)>=delta)&(2*abs(ex-ey)/(ex+ey)>=delta);out[f'rate{int(delta*100)}']=np.mean(good)
  rev.append(out)
  if s==103 and {a,b}=={'SAEPSO','SAC_SAPSO'}:
   xx=x[x.run==1].iloc[0];yy=y[y.run==1].iloc[0]
   for z in [xx,yy]:examples.append(dict(scenario=s,method=z.method,run=1,length=z.length,fR=z.fR,eta=z.fR/z.fL,harmonic=8*z.fL/z.fR))
pd.DataFrame(rev).to_csv(O/'reversal_materiality.csv',index=False);pd.DataFrame(examples).to_csv(O/'reversal_examples.csv',index=False)
# Rescore saved paths only. These rankings are not outcomes of reoptimization.
weights=[]
for s in maps:
 for exposure in ['accumulated','normalized']:
  for k in [0,1,2,3]:
   for scale in [.5,1,2]:
    q=A[A.scenario==s].copy();terms=q[['fL','fS','fH','fR']].to_numpy(copy=True)
    if exposure=='normalized':terms[:,3]/=terms[:,0]
    w=np.ones(4);w[k]=scale;q['reJ']=terms@w/w.sum();med=q.groupby('method').reJ.median().sort_values()
    weights.append(dict(scenario=s,exposure=exposure,term=k,scale=scale,lowestMedian=med.index[0],median=med.iloc[0],ordering='|'.join(med.index)))
pd.DataFrame(weights).to_csv(O/'rescore_weights.csv',index=False)
# Route bootstrap bands avoid seven heavy overlapping curves.
rng=np.random.default_rng(20261112);C=pd.DataFrame(curves);cr=[]
for (s,m),g in C.groupby(['scenario','method']):
 matrix=g.pivot(index='run',columns='c',values='F');v=matrix.to_numpy();boot=np.mean(v[rng.integers(0,30,(2000,30))],axis=1);lo,hi=np.quantile(boot,[.025,.975],axis=0)
 cr.extend([dict(scenario=s,method=m,c=c,mean=mean,low=l,high=h) for c,mean,l,h in zip(matrix.columns,v.mean(axis=0),lo,hi)])
pd.DataFrame(cr).to_csv(O/'cdf_intervals.csv',index=False)
summary=dict(originalPaths=len(A),minimumVertical=float(A.clearance.min()),endpoints=en,reversal={f'rate{k}':np.quantile(pd.DataFrame(rev)[f'rate{k}'],[0,.25,.5,.75,1]).tolist() for k in [0,1,5,10]},maps=[])
G=pd.DataFrame(rows)
for s in maps:
 a=A[A.scenario==s];g=G[G.scenario==s];summary['maps'].append(dict(scenario=s,medianVerticalSpatialRatio=[float(g.ratioLow.median()),float(g.ratioHigh.median())],trimMinimumMedians={str(k):float(g[f'trim{k}'].median()) for k in [1,5,10]},slopeRatioSpearman=float(spearmanr(g.slopeNorm,g.ratioLow).statistic),searchMedianSeconds=float(a.seconds.median()),assessmentMedianSeconds=float(a.seconds_assessment.median())))
(O/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
