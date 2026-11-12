"""Independently evaluate the 2,160 additional saved polylines; no search is launched."""
from pathlib import Path
from collections import defaultdict
import argparse,csv,json,hashlib
import numpy as np
from scipy.io import loadmat
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import fisher_exact,binomtest
from clearance_profiles import from_route
from distance_bounds import Terrain,distance
ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--euclidean',choices=['saved','examples','all'],default='saved');args=p.parse_args()
def holm(values):
 order=np.argsort(values);a=np.empty(len(values));prev=0.
 for rank,i in enumerate(order):prev=max(prev,min(1.,values[i]*(len(values)-rank)));a[i]=prev
 return a

sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for f,h in json.loads((ROOT/'SHA256.json').read_text()).items():assert sha(ROOT/f)==h,f
read=lambda p:list(csv.DictReader(p.open()))
envs={s:loadmat(ROOT/f'data/maps/scenario{s}.mat',simplify_cells=True)['env'] for s in [101,102,103]}
surfaces={s:RegularGridInterpolator((e['terrainY'][:,0],e['terrainX'][0]),e['terrainGrid']) for s,e in envs.items()};terrains={};totals={};recomputed=0
for cohort,expected in [('endpoint',1260),('constraint',360),('sensitivity',540)]:
 folder=ROOT/'additional'/cohort;rows=read(folder/'runs.csv');assert len(rows)==expected
 ledger=json.loads((folder/('protocol.json' if cohort=='sensitivity' else 'cases.json')).read_text())
 keyname='id' if cohort=='sensitivity' else ('caseId' if cohort=='endpoint' else 'scenario')
 cases={int(c[keyname]):c for c in ledger['cases']}
 methods=ledger['methods'] if cohort=='endpoint' else ['PSO_Standard','LSHADE']
 runs=range(1,11) if cohort=='sensitivity' else range(1,31)
 expectedkeys={(cid,m,c,run) for cid,case in cases.items() for m in methods for c in (['vertical','euclidean'] if cohort=='constraint' else [case.get('condition',case.get('family'))]) for run in runs}
 actualkeys={(int(r['case']),r['method'],r['condition'],int(r['run'])) for r in rows}
 assert len(rows)==len(actualkeys) and actualkeys==expectedkeys
 paths=defaultdict(list)
 for r in read(folder/'vertices.csv'):paths[r['id']].append((int(r['vertex']),[float(r[k]) for k in ['x','y','z']]))
 assert set(paths)=={r['id'] for r in rows};groups=defaultdict(list);nvertical=nbelow=nboth=nbinding=nunresolved=0
 for r in rows:
  case=cases[int(r['case'])];s=int(r['map']);assert s==int(case.get('map',case.get('mapId',case.get('scenario'))))
  P=np.array([p for _,p in sorted(paths[r['id']])]);e=dict(envs[s]);e.update(startPoint=P[0],goalPoint=P[-1]);h=float(r['required']);assert h==(case['h'] if cohort=='sensitivity' else 8)
  assert np.array_equal(P[0],case['start']) and np.array_equal(P[-1],case['goal'])
  assert int(r['waypoints'])==(case['n'] if cohort=='sensitivity' else 5)
  mi=methods.index(r['method'])+1;run=int(r['run'])
  seed=(5000000+(int(r['case'])-401)*10000+10*run+mi) if cohort=='sensitivity' else ((1000000+(int(r['case'])-201)*100000+10*run+mi) if cohort=='endpoint' else (3000000+(s-101)*100000+10*run+mi))
  assert int(r['seed'])==seed and int(r['searchFE'])==100100 and int(r['reportFE'])>=1
  assert P.shape==(int(r['waypoints'])+2,3) and np.isfinite(P).all()
  profile=from_route(P,e);assert abs(profile.minimum-float(r['verticalMinimum']))<1e-8
  assert abs(profile.length-float(r['length']))<1e-8
  D=np.diff(P,axis=0);ell=np.linalg.norm(D,axis=1)
  valid=np.all(P>=0) and np.all(P<=e['mapSize']) and np.all(ell>16*np.spacing(float(np.max(e['mapSize']))))
  vertical=bool(valid and profile.minimum>=h);assert vertical==bool(int(r['verticalFeasible']));nvertical+=vertical
  if vertical:
   theta=np.arctan2(np.linalg.norm(np.cross(D[:-1],D[1:]),axis=1),np.sum(D[:-1]*D[1:],axis=1))
   terms=np.array([profile.length/profile.direct_length,np.mean((theta/np.pi)**2),np.sum(ell*(P[:-1,2]+P[1:,2])/2)/profile.length/e['mapSize'][2],profile.inverse_profile(64)*profile.length/profile.direct_length])
   np.testing.assert_allclose(terms,[float(r[k]) for k in ['fL','fS','fH','fR']],rtol=0,atol=1e-7)
   assert abs(terms.mean()-float(r['J']))<1e-7
  rho=profile.minimum/h
  if cohort=='constraint' and r['condition']=='euclidean' and profile.minimum>0:
   search_lo,search_hi=float(r['searchDistanceLower']),float(r['searchDistanceUpper'])
   assert search_lo<=float(r['distanceUpper'])+1e-7 and float(r['distanceLower'])<=search_hi+1e-7
   rho=min(rho,search_lo/h)
  search_feasible=bool(valid and rho>=1)
  if cohort=='constraint':assert search_feasible==bool(int(r['searchFeasible']))
  violation=max(0,1-rho)
  score=(float(r['J'])/(1+float(r['J']))) if search_feasible else (1+violation/(1+violation) if valid else 2.)
  assert abs(score-float(r['score']))<1e-7
  if cohort=='endpoint':
   ends=P[[0,-1],2]-surfaces[s](P[[0,-1]][:,[1,0]]);binding=bool(abs(profile.minimum-min(ends))<=1e-8)
   assert binding==bool(int(r['endpointBinding']));nbinding+=binding
  else:
   lo,hi=float(r['distanceLower']),float(r['distanceUpper']);assert 0<=lo<=hi and hi-lo<=.00100001
   below=hi<h;both=vertical and lo>=h;nbelow+=below;nboth+=both;nunresolved+=lo<h<=hi
   q=np.array([float(r[k]) for k in ['pathX','pathY','pathZ']]);t=np.array([float(r[k]) for k in ['surfaceX','surfaceY','surfaceZ']])
   u=np.clip(np.sum((q-P[:-1])*D,axis=1)/np.sum(D*D,axis=1),0,1)
   assert np.min(np.linalg.norm(q-(P[:-1]+u[:,None]*D),axis=1))<1e-8
   assert abs(float(surfaces[s](t[[1,0]])[0])-t[2])<1e-7
   assert abs(np.linalg.norm(q-t)+1e-9-float(r['nativeUpper']))<1e-8
   assert float(r['nativeLower'])<=hi+1e-7 and lo<=float(r['nativeUpper'])+1e-7
   if args.euclidean=='all' or (args.euclidean=='examples' and int(r['run'])==1):
    if s not in terrains:terrains[s]=Terrain(e)
    out=distance(P,terrains[s],tolerance=1e-4)
    if out['lower']<h<=out['upper']:out=distance(P,terrains[s],tolerance=1e-7)
    assert out['converged'] and out['lower']<=hi+1e-7 and lo<=out['upper']+1e-7;recomputed+=1
  groups[r['case'],r['method']].append(r)
 totals[cohort]=dict(routes=len(rows),verticallyFeasible=nvertical,belowSpatial=nbelow if cohort!='endpoint' else None,bothMargins=nboth if cohort!='endpoint' else None,distanceUnresolved=nunresolved if cohort!='endpoint' else None,endpointBinding=nbinding if cohort=='endpoint' else None)
 if cohort=='endpoint':
  tests=[];summary=read(folder/'original_assessment/summary.csv');saved=read(folder/'original_assessment/endpoint_binding_tests.csv')
  for row in summary:
   g=groups[row['case'],row['method']];good=[r for r in g if int(r['verticalFeasible'])]
   assert len(g)==int(row['runs']) and len(good)==int(row['feasible']) and sum(int(r['endpointBinding']) for r in good)==int(row['boundAttained'])
   if good:np.testing.assert_allclose(np.quantile([float(r['J']) for r in good],[.25,.5,.75]),[float(row[k]) for k in ['J_q1','J_median','J_q3']],rtol=0,atol=1e-12)
  for row in saved:
   bycond={c:[r for r in rows if r['map']==row['map'] and r['method']==row['method'] and r['condition']==c and int(r['verticalFeasible'])] for c in ['near','high']}
   counts=[[sum(int(r['endpointBinding']) for r in bycond[c]),sum(1-int(r['endpointBinding']) for r in bycond[c])] for c in ['near','high']]
   tests.append(float(fisher_exact(counts,alternative='two-sided').pvalue))
  np.testing.assert_allclose(tests,[float(r['p']) for r in saved],rtol=1e-10,atol=1e-14)
  np.testing.assert_allclose(holm(tests),[float(r['holm_p']) for r in saved],rtol=1e-10,atol=1e-14)
 if cohort=='constraint':
  initial=np.load(folder/'initial_populations.npz');index={(r['map'],r['method'],r['condition'],r['run']):r for r in rows}
  native_interval_checks=[]
  for s in ['101','102','103']:
   for m in ['PSO_Standard','LSHADE']:
    four=np.zeros((2,2),int)
    coarse_unresolved=search_confirmed=independent_unresolved=0
    for run in map(str,range(1,31)):
     a,b=[index[s,m,c,run] for c in ['vertical','euclidean']]
     np.testing.assert_array_equal(initial[a['id']],initial[b['id']])
     aok=int(a['verticalFeasible']) and float(a['distanceLower'])>=8;bok=int(b['verticalFeasible']) and float(b['distanceLower'])>=8
     four[int(aok),int(bok)]+=1
     coarse_unresolved+=float(b['nativeLower'])<8<=float(b['nativeUpper'])
     search_confirmed+=int(b['searchFeasible'])
     independent_unresolved+=float(b['distanceLower'])<8<=float(b['distanceUpper'])
    assert four.sum()==30
    native_interval_checks.append(dict(map=int(s),method=m,nativePostHocUnresolved=coarse_unresolved,nativeSearchConfirmed=search_confirmed,independentUnresolved=independent_unresolved))
  totals[cohort]['euclideanArmIntervals']=dict(nativePostHocTolerance=1e-5,thresholdRefinementTolerance=1e-7,conditions=native_interval_checks)
  # Recompute all six test families and pointwise paired bootstrap summaries.
  from paired_statistics import analyze
  converted=[]
  integer={'map','case','run','seed','waypoints','verticalFeasible','searchFeasible','geometryCalls','thresholdAmbiguousCalls','searchFE','reportFE'}
  strings={'id','method','condition','rawSHA256','protocolHash'}
  for row in rows:
   r={k:(v if k in strings else (int(v) if k in integer else float(v))) for k,v in row.items()}
   r['scenario']=r['map'];r['bothMargins']=int(r['verticalFeasible'] and r['distanceLower']>=8);r['distanceUnresolved']=int(r['distanceLower']<8<=r['distanceUpper']);converted.append(r)
  dest=ROOT/'recomputed_paired_statistics';dest.mkdir(exist_ok=True);analyze(converted,dest)
  observed=read(dest/'paired_summary.csv');saved=read(folder/'original_assessment/paired_summary.csv');assert len(observed)==len(saved)==6
  for a,b in zip(observed,saved):
   assert a['scenario']==b['scenario'] and a['method']==b['method']
   for k in a:
    if k in ['scenario','method']:continue
    np.testing.assert_allclose(float(a[k]),float(b[k]),rtol=1e-10,atol=1e-10,equal_nan=True)
  boundary=json.loads((ROOT/'revised_analysis/constraint_boundary_summary.json').read_text())
  assert boundary['verified'] and len(boundary['conditions'])==6
  for result in boundary['conditions']:
   g=[r for r in rows if int(r['map'])==result['scenario'] and r['method']==result['method'] and r['condition']=='euclidean']
   assert len(g)==result['runs']==30
   lo=np.array([float(r['distanceLower']) for r in g]);hi=np.array([float(r['distanceUpper']) for r in g]);amb=np.array([int(r['thresholdAmbiguousCalls']) for r in g])
   assert int((hi<=8.001).sum())==result['upperWithin001'] and int((hi<=8.002).sum())==result['upperWithin002']
   assert min(lo)-8==result['minimumLowerExcess'] and max(hi)-8==result['maximumUpperExcess']
   assert np.median(amb)==result['ambiguousEvaluationMedian'] and min(amb)==result['ambiguousEvaluationMinimum'] and max(amb)==result['ambiguousEvaluationMaximum']
 if cohort=='sensitivity':
  # Endpoint screening is independent of the saved optimizer outcomes. The
  # short segment's lower bound bounds the endpoint; adding its length to its
  # upper bound gives a conservative endpoint upper bound (1-Lipschitz).
  endpoint_records=json.loads((ROOT/'revised_analysis/sensitivity_endpoint_distances.json').read_text())
  endpoint_by_key={(int(r['case']),r['endpoint']):r for r in endpoint_records}
  assert len(endpoint_by_key)==54
  endpoint_cache={};endpoint_failed=set()
  for cid,case in cases.items():
   s=int(case['map'])
   if s not in terrains:terrains[s]=Terrain(envs[s])
   for which in ['start','goal']:
    point=np.array(case[which],dtype=float);key=(s,*point)
    if key not in endpoint_cache:
     eps=1e-7;out=distance(np.vstack([point,point+[0,0,eps]]),terrains[s],tolerance=1e-6)
     assert out['converged'];endpoint_cache[key]=(out['lower'],out['upper']+eps)
    lo,hi=endpoint_cache[key];saved=endpoint_by_key[cid,which]
    assert saved['coordinates']==case[which] and saved['required']==case['h']
    assert lo<=saved['upper']+1e-7 and saved['lower']<=hi+1e-7
    assert hi<case['h'] or lo>=case['h']
    assert (hi<case['h'])==saved['below']
    if hi<case['h']:endpoint_failed.add(cid)
  assert endpoint_failed=={407,408,425,426,427}
  witnesses=json.loads((ROOT/'revised_analysis/sensitivity_feasible_witnesses.json').read_text())
  assert witnesses['verified'] and {r['case'] for r in witnesses['paths']}==set(cases)-endpoint_failed
  for witness in witnesses['paths']:
   case=cases[witness['case']];s=int(case['map']);P=np.array(witness['vertices'])
   assert P.shape==(case['n']+2,3) and np.array_equal(P[0],case['start']) and np.array_equal(P[-1],case['goal'])
   assert np.isfinite(P).all() and np.all(P>=0) and np.all(P<=envs[s]['mapSize'])
   assert np.all(np.linalg.norm(np.diff(P,axis=0),axis=1)>16*np.spacing(float(np.max(envs[s]['mapSize']))))
   e=dict(envs[s]);e.update(startPoint=P[0],goalPoint=P[-1]);profile=from_route(P,e)
   assert profile.minimum>=case['h'] and abs(profile.minimum-witness['verticalMinimum'])<1e-8
   out=distance(P,terrains[s],tolerance=1e-5)
   assert out['converged'] and out['lower']>=case['h']
   assert out['lower']<=witness['distanceUpper']+1e-7 and witness['distanceLower']<=out['upper']+1e-7
  totals[cohort]['verifiedExistenceWitnesses']=22
  forced=[r for r in rows if int(r['case']) in endpoint_failed]
  assert len(forced)==100 and all(float(r['distanceUpper'])<float(r['required']) for r in forced)
  totals[cohort]['endpointExcludedConfigurations']=sorted(endpoint_failed)
  totals[cohort]['unavoidableEndpointBreaches']=len(forced)
  eligible=[r for r in rows if int(r['case']) not in endpoint_failed]
  totals[cohort]['endpointAdmissibleRoutes']=len(eligible)
  totals[cohort]['belowSpatialWithAdmissibleEndpoints']=sum(float(r['distanceUpper'])<float(r['required']) for r in eligible)
  summary=read(ROOT/'revised_analysis/sensitivity_summary.csv');assert len(summary)==54
  for row in summary:
   g=groups[row['case'],row['method']];good=[r for r in g if int(r['verticalFeasible'])]
   assert len(g)==int(row['runs'])==10 and len(good)==int(row['verticalFeasible'])
   assert (row['endpointThresholdMet']=='True')==(int(row['case']) not in endpoint_failed)
   assert sum(float(r['distanceUpper'])<float(r['required']) for r in g)==int(row['belowSpatial'])
   assert sum(int(r['verticalFeasible']) and float(r['distanceLower'])>=float(r['required']) for r in g)==int(row['bothMargins'])
   assert sum(float(r['distanceLower'])<float(r['required'])<=float(r['distanceUpper']) for r in g)==int(row['unresolved'])
   if good:np.testing.assert_allclose(np.quantile([float(r['J']) for r in good],[.25,.5,.75]),[float(row[k]) for k in ['JQ1','JMedian','JQ3']],rtol=0,atol=1e-12)
   for source,target in [('distanceLower','dMedianLow'),('distanceUpper','dMedianHigh'),('seconds','secondsMedian')]:assert abs(np.median([float(r[source]) for r in g])-float(row[target]))<1e-9
   def wilson(k,n):
    z=1.959963984540054;c=(k/n+z*z/(2*n))/(1+z*z/n);h=z*np.sqrt((k/n)*(1-k/n)/n+z*z/(4*n*n))/(1+z*z/n);return max(0,c-h),min(1,c+h)
   lo=wilson(int(row['belowSpatial']),10)[0];hi=wilson(int(row['belowSpatial'])+int(row['unresolved']),10)[1]
   assert abs(lo-float(row['belowWilsonLow']))<1e-12 and abs(hi-float(row['belowWilsonHigh']))<1e-12
 print(cohort,totals[cohort],flush=True)
assert sum(x['routes'] for x in totals.values())==2160
report=dict(verified=True,additional=totals,independentlyRecomputedEuclideanDistances=recomputed,initialPopulationPairs=180,endpointFisherHolmContrasts=21,pairedMcNemarHolmContrasts=6,pairedBootstrapSummaries=True,sensitivityConditions=54)
(ROOT/'additional_reanalysis_report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
