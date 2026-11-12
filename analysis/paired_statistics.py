import csv,math
import numpy as np
from scipy.stats import binomtest,beta
METHODS=['PSO_Standard','LSHADE']
def writecsv(path,rows):
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def holm(p):
 order=np.argsort(p);out=np.zeros(len(p));prev=0
 for rank,i in enumerate(order):prev=max(prev,min(1.,p[i]*(len(p)-rank)));out[i]=prev
 return out

def paired_rate_interval(n01,n10,n):
 # Each discordance probability has a97.5% exact binomial interval.
 # Bonferroni gives at least95% joint coverage; differencing gives a
 # conservative interval for p01-p10, including all-zero/all-one samples.
 def cp(k):
  return (0. if k==0 else float(beta.ppf(.0125,k,n-k+1)),1. if k==n else float(beta.ppf(.9875,k+1,n-k)))
 l01,u01=cp(n01);l10,u10=cp(n10)
 return l01-u10,u01-l10

def analyze(rows,folder):
 rng=np.random.default_rng(20261112);geometry_rng=np.random.default_rng(2026111208);summaries=[];look={(r['scenario'],r['method'],r['condition'],r['run']):r for r in rows}
 for s in [101,102,103]:
  for method in METHODS:
   pairs=[(look[s,method,'vertical',r],look[s,method,'euclidean',r]) for r in range(1,31)]
   counts=[sum(a['bothMargins']==i and b['bothMargins']==j for a,b in pairs) for i,j in [(0,0),(0,1),(1,0),(1,1)]]
   n00,n01,n10,n11=counts;discord=n01+n10;p=float(binomtest(n01,discord,.5).pvalue) if discord else 1.
   out=dict(scenario=s,method=method,pairs=30,bothFail=n00,euclideanOnly=n01,verticalOnly=n10,bothSucceed=n11,verticalRate=(n10+n11)/30,euclideanRate=(n01+n11)/30,mcnemarP=p)
   for metric in ['J','fL','fS','fH','fR','length','seconds']:
    vals=np.array([b[metric]-a[metric] for a,b in pairs if np.isfinite(a[metric]) and np.isfinite(b[metric])]);out[metric+'Pairs']=len(vals)
    if len(vals):
     boot=np.median(vals[rng.integers(0,len(vals),(10000,len(vals)))],axis=1);out[metric+'DeltaMedian']=float(np.median(vals));out[metric+'Low'],out[metric+'High']=map(float,np.quantile(boot,[.025,.975]))
    else:out[metric+'DeltaMedian']=out[metric+'Low']=out[metric+'High']=math.nan
   out['verticalUnresolved']=sum(a['distanceUnresolved'] for a,b in pairs)
   out['euclideanUnresolved']=sum(b['distanceUnresolved'] for a,b in pairs)
   out['verticalGeometryCallsMedian']=float(np.median([a['geometryCalls'] for a,b in pairs]))
   out['euclideanGeometryCallsMedian']=float(np.median([b['geometryCalls'] for a,b in pairs]))
   out['euclideanAmbiguousCallsMedian']=float(np.median([b['thresholdAmbiguousCalls'] for a,b in pairs]))
   dlo=np.array([b['distanceLower']-a['distanceUpper'] for a,b in pairs]);dhi=np.array([b['distanceUpper']-a['distanceLower'] for a,b in pairs])
   ix=geometry_rng.integers(0,30,(10000,30))
   out['distanceDeltaMedianLower']=float(np.median(dlo));out['distanceDeltaMedianUpper']=float(np.median(dhi))
   out['distanceDeltaCILower']=float(np.quantile(np.median(dlo[ix],axis=1),.025));out['distanceDeltaCIUpper']=float(np.quantile(np.median(dhi[ix],axis=1),.975))
   success=np.array([b['bothMargins']-a['bothMargins'] for a,b in pairs]);out['successRateDifference']=float(success.mean())
   out['successRateDifferenceLow'],out['successRateDifferenceHigh']=paired_rate_interval(n01,n10,30)
   summaries.append(out)
 for r,p in zip(summaries,holm([r['mcnemarP'] for r in summaries])):r['holmP']=float(p)
 writecsv(folder/'paired_summary.csv',summaries)
 lines=['# Verified clearance-definition intervention','', 'The360 primary routes form180 seed-matched pairs on the three original maps and endpoints. The objective and search budget are unchanged. Initial candidate populations match within every pair. Primary inference concerns both modeled margins, with six exact McNemar contrasts corrected by Holm. Secondary intervals are pointwise route-pair bootstrap estimates. Cost estimates condition on pairs with a defined original objective in both arms.','', '| Map | Method | Vertical-arm success /30 | Euclidean-arm success /30 | Holm p | Median deltaJ | Cost pairs |','|---|---|---:|---:|---:|---:|---:|']
 for r in summaries:lines.append(f"| {r['scenario']} | {r['method']} | {r['verticalOnly']+r['bothSucceed']} | {r['euclideanOnly']+r['bothSucceed']} | {r['holmP']:.5g} | {r['JDeltaMedian']:.6g} | {r['JPairs']} |")
 lines+=['','Success refers to distance from the continuous bilinear height field together with the original vertical requirement. It is not physical flight validation. These results are conditional on two methods and three fixed missions. Original630-route and separate endpoint cohorts remain distinct.']
 (folder/'RESULTS.md').write_text('\n'.join(lines)+'\n')
