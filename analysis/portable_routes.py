from pathlib import Path
import csv
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parent
VERTICES={}
for r in csv.DictReader((ROOT/'data/returned_routes.csv').open()):
 k=(int(r['scenario']),r['method'],int(r['run']));VERTICES.setdefault(k,[]).append((int(r['vertex']),[float(r[n]) for n in ['x','y','z']]))
ENVS={s:loadmat(ROOT/f'data/maps/scenario{s}.mat',simplify_cells=True)['env'] for s in [101,102,103]}
def read_original_route(s,m,r):
 return {'path':np.array([p for _,p in sorted(VERTICES[int(s),m,int(r)])]),'env':ENVS[int(s)]}
