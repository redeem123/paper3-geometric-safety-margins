"""Secondary effect-size estimates from the unchanged 630 independent runs."""
from pathlib import Path
import csv,hashlib,json
import numpy as np
from scipy.stats import mannwhitneyu

root=Path(__file__).resolve().parent
source=root/'verified_runs.csv'
rows=list(csv.DictReader(source.open(newline='')))
assert len(rows)==630 and all(int(r['band'])==3 for r in rows)
comparisons={(int(r['scenario']),r['comparator']):r for r in csv.DictReader((root/'comparisons.csv').open())}
methods=[('PSO_Standard','PSO'),('SAEPSO','SAEPSO'),('SAC_SAPSO','SAC-SAPSO'),('RLAM_PSO','RLAM-PSO'),('LSHADE','L-SHADE'),('GLEET','GLEET')]
scenarios=[(101,'Urban'),(102,'Suburban'),(103,'Mountain')]
B=20000
streams=np.random.SeedSequence(20261111).spawn(18)
result=[]
for s,name in scenarios:
    x=np.array([float(r['J']) for r in rows if int(r['scenario'])==s and r['method']=='AFSACPSO'])
    assert len(x)==30
    for method,display in methods:
        y=np.array([float(r['J']) for r in rows if int(r['scenario'])==s and r['method']==method])
        assert len(y)==30
        a=float(np.mean((x[:,None]<y[None,:])+.5*(x[:,None]==y[None,:])))
        u=mannwhitneyu(x,y,alternative='two-sided',method='asymptotic').statistic
        assert abs(a-(1-u/(len(x)*len(y))))<1e-12
        rng=np.random.default_rng(streams[len(result)])
        differences=[]
        for start in range(0,B,500):
            n=min(500,B-start)
            xb=x[rng.integers(0,len(x),size=(n,len(x)))];yb=y[rng.integers(0,len(y),size=(n,len(y)))]
            differences.append(np.median(xb,axis=1)-np.median(yb,axis=1))
        dci=np.quantile(np.concatenate(differences),[.025,.975])
        result.append(dict(scenario=s,map=name,method=method,display=display,nAF=len(x),nComparator=len(y),
            deltaMedian=float(np.median(x)-np.median(y)),deltaLow=float(dci[0]),deltaHigh=float(dci[1]),
            probabilityAF=a,
            holmScore=float(comparisons[(s,method)]['holm_score'])))
with (root/'secondary_effect_sizes.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(result[0]));writer.writeheader();writer.writerows(result)
metadata=dict(analysis='Secondary descriptive analysis, added after the primary study',sourceSHA256=hashlib.sha256(source.read_bytes()).hexdigest(),bootstrapReplicates=B,bootstrapMethod='Independent within-method resampling of 30 runs; percentile pointwise 95% intervals',randomGenerator='NumPy default_rng from SeedSequence(20261111), one spawned stream per map/comparator',deltaDefinition='median(J_AFSACPSO) - median(J_comparator)',probabilityDefinition='P(J_AFSACPSO < J_comparator) + 0.5 P(tie)',multiplicity='Intervals are pointwise and not adjusted for simultaneous coverage; original Holm-adjusted score tests are unchanged',results=result)
(root/'secondary_effect_sizes.json').write_text(json.dumps(metadata,indent=2)+'\n')
lines=[r'\begin{table}[htbp]',r'\centering\small\setlength{\tabcolsep}{4pt}',
       r'\caption{Secondary effect-size estimates for AFSACPSO versus each comparator, using 30 independent runs per method and map. $\Delta_J$ entries give median differences (AFSACPSO minus comparator) [pointwise 95\% bootstrap intervals]. $\widehat A$ gives descriptive estimates of the probability of a lower AFSACPSO cost with half credit for ties. The final column retains the original Holm-adjusted score-test $p$ values.}',
       r'\label{tab:effect-sizes}',r'\begin{tabular}{llccc}\toprule',r'Map & Comparator & $\Delta_J$ & $\widehat A$ & Adjusted $p$\\\midrule']
for j,r in enumerate(result):
    if j and j%6==0:lines.append(r'\midrule')
    name=r['map'] if j%6==0 else ''
    delta=f"{r['deltaMedian']:.5f} [{r['deltaLow']:.5f}, {r['deltaHigh']:.5f}]"
    prob=f"{r['probabilityAF']:.3f}"
    pv=r['holmScore'];ex=int(np.floor(np.log10(pv)));coeff=pv/10**ex
    ps=f'${coeff:.2f}\\times10^{{{ex}}}$' if pv<.001 else f'{pv:.4f}'
    lines.append(f"{name} & {r['display']} & {delta} & {prob} & {ps}"+r'\\')
lines += [r'\bottomrule\end{tabular}',r'\end{table}','']
(root/'tables/effect_sizes.tex').write_text('\n'.join(lines))
for r in result:
    if r['method'] in ('SAC_SAPSO','LSHADE'):print(r['map'],r['display'],'Delta',f"{r['deltaMedian']:.7f}",'CI',round(r['deltaLow'],7),round(r['deltaHigh'],7),'A',round(r['probabilityAF'],4))
print('Wrote all 18 secondary contrasts. Original result and primary-test files unchanged.')
