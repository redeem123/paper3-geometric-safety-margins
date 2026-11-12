"""Arc-length clearance distribution from bilinear-cell quadratics.

The cellwise terrain representation follows the previously used continuous
minimum analysis. These routines derive profile summaries, not a new optimizer.
"""
from dataclasses import dataclass
import numpy as np
from scipy.interpolate import RegularGridInterpolator


def sublevel_measure(a,b,c,thresholds):
    """Lebesgue measure in u in [0,1] where a*u*u+b*u+c <= threshold."""
    a,b,c=np.broadcast_arrays(np.asarray(a,float),np.asarray(b,float),np.asarray(c,float))
    a=a.ravel();b=b.ravel();c=c.ravel();t=np.atleast_1d(np.asarray(thresholds,float))
    result=np.zeros((len(a),len(t)))
    # Near-flat intervals are not atoms: arbitrarily small nonzero coefficients
    # can still give an O(1) CDF change at thresholds inside their support.
    quadratic=a!=0
    linear=(~quadratic)&(b!=0)
    constant=~(quadratic|linear)
    if np.any(constant):result[constant]=(c[constant,None]<=t[None,:])
    if np.any(linear):
        cut=(t[None,:]-c[linear,None])/b[linear,None]
        fraction=np.clip(cut,0,1)
        result[linear]=np.where(b[linear,None]>0,fraction,1-fraction)
    if np.any(quadratic):
        A=a[quadratic,None];B=b[quadratic,None];C=c[quadratic,None]-t[None,:]
        disc=B*B-4*A*C;real=disc>0
        sqrt=np.sqrt(np.maximum(disc,0))
        q=-.5*(B+np.copysign(sqrt,B))
        r1=np.zeros_like(C);r2=np.zeros_like(C)
        np.divide(q,A,out=r1,where=real)
        np.divide(C,q,out=r2,where=real)
        low=np.minimum(r1,r2);high=np.maximum(r1,r2)
        inside=np.maximum(0,np.minimum(1,high)-np.maximum(0,low))
        values=np.where(A>0,inside,1-inside)
        # A positive double-root parabola contributes only a measure-zero point;
        # negative curvature contributes the full interval at/above its maximum.
        result[quadratic]=np.where(real,values,(A<0).astype(float))
    return np.clip(result,0,1)


@dataclass
class ClearanceProfile:
    a:np.ndarray
    b:np.ndarray
    c:np.ndarray
    weights:np.ndarray  # ell_segment * dt / L
    length:float
    direct_length:float
    minimum:float
    maximum:float

    def cdf(self,thresholds):
        out=self.weights@sublevel_measure(self.a,self.b,self.c,thresholds)
        assert np.min(out)>=-1e-12 and np.max(out)<=1+1e-12
        return np.clip(out,0,1)

    def quantiles(self,probabilities):
        p=np.atleast_1d(np.asarray(probabilities,float));assert np.all((p>0)&(p<1))
        low=np.full(p.shape,self.minimum);high=np.full(p.shape,self.maximum)
        for _ in range(44):
            mid=(low+high)/2;left=self.cdf(mid)>=p
            high=np.where(left,mid,high);low=np.where(left,low,mid)
        return high

    def inverse_profile(self,order=128,required=8.):
        assert self.minimum>0 and required>0
        nodes,w=np.polynomial.legendre.leggauss(order);u=(nodes+1)/2
        h=self.a[:,None]*u**2+self.b[:,None]*u+self.c[:,None]
        assert np.min(h)>0
        return float(required*np.sum(self.weights*((1/h)@w)/2))


def from_route(path,env):
    path=np.asarray(path,float);assert path.ndim==2 and path.shape[1]==3 and np.isfinite(path).all()
    terrain=RegularGridInterpolator((env['terrainY'][:,0],env['terrainX'][0,:]),env['terrainGrid'])
    d=np.diff(path,axis=0);ell=np.linalg.norm(d,axis=1);L=float(ell.sum());L0=float(np.linalg.norm(path[-1]-path[0]))
    assert np.all(ell>0) and L0>0
    aa=[];bb=[];cc=[];ww=[]
    for q,v,length in zip(path[:-1],d,ell):
        cuts=[np.array([0.,1.])]
        for dim,axis in ((0,terrain.grid[1]),(1,terrain.grid[0])):
            if v[dim]!=0:
                t=(axis-q[dim])/v[dim];cuts.append(t[(t>0)&(t<1)])
        t=np.unique(np.concatenate(cuts));left=t[:-1];width=np.diff(t)
        def h(at):
            p=q+at[:,None]*v
            return p[:,2]-terrain(p[:,[1,0]])
        c0,c1,cm=h(left),h(left+width),h(left+width/2)
        A=2*(c0+c1-2*cm);B=c1-c0-A
        aa.extend(A);bb.extend(B);cc.extend(c0);ww.extend(length*width/L)
    a=np.asarray(aa);b=np.asarray(bb);c=np.asarray(cc);weights=np.asarray(ww)
    assert abs(weights.sum()-1)<1e-12
    ends=np.stack([c,a+b+c]);mins=np.min(ends,axis=0);maxs=np.max(ends,axis=0)
    u=np.zeros_like(a);np.divide(-b,2*a,out=u,where=a!=0)
    inside=(u>0)&(u<1)&(a!=0);vertex=a*u*u+b*u+c
    mins=np.where(inside&(a>0),np.minimum(mins,vertex),mins)
    maxs=np.where(inside&(a<0),np.maximum(maxs,vertex),maxs)
    return ClearanceProfile(a,b,c,weights,L,L0,float(mins.min()),float(maxs.max()))
