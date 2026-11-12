"""Geometric branch-and-bound for polyline versus a bilinear height field.

AABB bounds enclose every bilinear patch. Point witnesses supply upper bounds.
Floating-point results include an explicit numerical allowance; this is not
interval-arithmetic certification. No optimizer or benchmark is launched.
"""
import heapq
import numpy as np
from triangle_bound import patch_triangle_bound
from scipy.interpolate import RegularGridInterpolator

GUARD=1e-9

def box_bounds(path,lo,hi):
    """Exact piecewise-quadratic segment/AABB distance, up to roundoff.

    Return a lower bound and the closest box point for every box, taking the
    minimum over all segments. Segment breakpoints are the six box-plane hits.
    """
    p=np.asarray(path,float);lo=np.atleast_2d(lo);hi=np.atleast_2d(hi)
    a=p[:-1];v=np.diff(p,axis=0);B=len(lo);N=len(a)
    hits=[]
    for bound in [lo,hi]:
        h=np.zeros((B,N,3));np.divide(bound[:,None,:]-a[None,:,:],v[None,:,:],out=h,where=v[None,:,:]!=0)
        hits.append(np.clip(h,0,1))
    t=np.sort(np.concatenate([np.zeros((B,N,1)),np.ones((B,N,1)),*hits],axis=2),axis=2)
    left=t[:,:,:-1];right=t[:,:,1:];mid=(left+right)/2
    q=a[None,:,None,:]+mid[:,:,:,None]*v[None,:,None,:]
    low=q<lo[:,None,None,:];high=q>hi[:,None,None,:];outside=low|high
    A=np.where(low,a[None,:,None,:]-lo[:,None,None,:],np.where(high,a[None,:,None,:]-hi[:,None,None,:],0))
    V=np.where(outside,v[None,:,None,:],0)
    den=np.sum(V*V,axis=3);num=-np.sum(A*V,axis=3)
    stationary=np.array(mid);np.divide(num,den,out=stationary,where=den>0)
    best=np.clip(stationary,left,right)
    q=a[None,:,None,:]+best[:,:,:,None]*v[None,:,None,:]
    boxq=np.minimum(hi[:,None,None,:],np.maximum(lo[:,None,None,:],q))
    dsq=np.sum((q-boxq)**2,axis=3)
    ix=np.argmin(dsq.reshape(B,-1),axis=1);seg=ix//7;interval=ix%7
    rows=np.arange(B);d=np.sqrt(dsq[rows,seg,interval])
    witness=boxq[rows,seg,interval]
    return np.maximum(0,d-GUARD),witness

def point_upper(path,points):
    a=path[:-1];v=np.diff(path,axis=0);points=np.atleast_2d(points)
    t=np.sum((points[:,None,:]-a)*v,axis=2)/np.sum(v*v,axis=1)
    t=np.clip(t,0,1);q=a+t[:,:,None]*v
    d=np.linalg.norm(points[:,None,:]-q,axis=2);i,j=np.unravel_index(np.argmin(d),d.shape)
    return float(d[i,j]+GUARD),q[i,j],points[i],int(j),float(t[i,j])

class Terrain:
    def __init__(self,env):
        self.x=np.asarray(env['terrainX'][0],float);self.y=np.asarray(env['terrainY'][:,0],float);self.H=np.asarray(env['terrainGrid'],float)
        self.interp=RegularGridInterpolator((self.y,self.x),self.H)
        self.nodes=[]
        def build(i0,i1,j0,j1):
            heights=self.H[i0:i1+1,j0:j1+1]
            index=len(self.nodes);self.nodes.append(None)
            lo=np.array([self.x[j0],self.y[i0],heights.min()]);hi=np.array([self.x[j1],self.y[i1],heights.max()])
            children=[]
            if i1-i0>1 or j1-j0>1:
                if (i1-i0)>=(j1-j0) and i1-i0>1:
                    c=(i0+i1)//2;children=[build(i0,c,j0,j1),build(c,i1,j0,j1)]
                else:
                    c=(j0+j1)//2;children=[build(i0,i1,j0,c),build(i0,i1,c,j1)]
            self.nodes[index]=(lo,hi,children,(i0,i1,j0,j1))
            return index
        build(0,len(self.y)-1,0,len(self.x)-1)
    def z(self,xy):
        xy=np.atleast_2d(xy);return self.interp(xy[:,[1,0]])
    def patch(self,node):
        i0,i1,j0,j1=self.nodes[node][3]
        return (self.x[j0],self.x[j1],self.y[i0],self.y[i1],self.H[i0,j0],self.H[i0,j1],self.H[i1,j0],self.H[i1,j1])

def patch_box(p):
    x0,x1,y0,y1,*h=p
    return np.array([x0,y0,min(h)]),np.array([x1,y1,max(h)])

def patch_value(p,xy):
    x0,x1,y0,y1,z00,z10,z01,z11=p;xy=np.atleast_2d(xy)
    u=np.clip((xy[:,0]-x0)/(x1-x0),0,1);v=np.clip((xy[:,1]-y0)/(y1-y0),0,1)
    return (1-u)*(1-v)*z00+u*(1-v)*z10+(1-u)*v*z01+u*v*z11

def split_patch(p):
    x0,x1,y0,y1,z00,z10,z01,z11=p;xm=(x0+x1)/2;ym=(y0+y1)/2
    zb=(z00+z10)/2;zt=(z01+z11)/2;zl=(z00+z01)/2;zr=(z10+z11)/2;zc=(z00+z10+z01+z11)/4
    return [(x0,xm,y0,ym,z00,zb,zl,zc),(xm,x1,y0,ym,zb,z10,zc,zr),
            (x0,xm,ym,y1,zl,zc,z01,zt),(xm,x1,ym,y1,zc,zr,zt,z11)]

def distance(path,terrain,tolerance=1e-3,max_pops=300000,initial_points=None):
    path=np.asarray(path,float);assert np.all(np.linalg.norm(np.diff(path,axis=0),axis=1)>0)
    # Endpoints and supplied exact surface-node witnesses initialize a valid UB.
    xy=path[:,0:2];points=np.column_stack([xy,terrain.z(xy)])
    if initial_points is not None:points=np.vstack([points,np.atleast_2d(initial_points)])
    upper,pw,sw,segment,t=point_upper(path,points)
    heap=[];counter=0;visited=0
    lo,hi,*_=terrain.nodes[0];lb,_=box_bounds(path,lo,hi)
    heapq.heappush(heap,(float(lb[0]),counter,0,0))
    while heap and upper-heap[0][0]>tolerance and visited<max_pops:
        lb,_,kind,obj=heapq.heappop(heap);visited+=1
        if lb>=upper:continue
        if kind==0:
            lo,hi,children,_=terrain.nodes[obj]
            if children:
                los=np.array([terrain.nodes[n][0] for n in children]);his=np.array([terrain.nodes[n][1] for n in children])
                lbs,boxq=box_bounds(path,los,his)
                sample=np.column_stack([boxq[:,:2],terrain.z(boxq[:,:2])]);u,q,r,j,tt=point_upper(path,sample)
                if u<upper:upper,pw,sw,segment,t=u,q,r,j,tt
                for n,l in zip(children,lbs):
                    if l<upper:counter+=1;heapq.heappush(heap,(float(l),counter,0,n))
                continue
            p=terrain.patch(obj)
        else:p=obj
        patches=split_patch(p);bounds=[patch_box(q) for q in patches]
        lbs,boxq=box_bounds(path,np.array([b[0] for b in bounds]),np.array([b[1] for b in bounds]))
        samples=[]
        for index,(pp,q) in enumerate(zip(patches,boxq)):
            triLower,triPoint=patch_triangle_bound(path,pp)
            lbs[index]=max(lbs[index],triLower)
            samples.append(np.column_stack([triPoint[None,:2],patch_value(pp,triPoint[None,:2])]))
            x0,x1,y0,y1,*_=pp
            xy=np.array([q[:2],[(x0+x1)/2,(y0+y1)/2],[x0,y0],[x1,y0],[x0,y1],[x1,y1]])
            samples.append(np.column_stack([xy,patch_value(pp,xy)]))
        u,q,r,j,tt=point_upper(path,np.concatenate(samples))
        if u<upper:upper,pw,sw,segment,t=u,q,r,j,tt
        for pp,l in zip(patches,lbs):
            if l<upper:counter+=1;heapq.heappush(heap,(float(l),counter,1,pp))
    lower=min(upper,heap[0][0]) if heap else upper
    assert lower<=upper+1e-10
    return dict(lower=float(lower),upper=float(upper),gap=float(upper-lower),converged=bool(upper-lower<=tolerance),
                pops=visited,pathPoint=pw.tolist(),surfacePoint=sw.tolist(),segment=segment,segmentFraction=t,tolerance=tolerance,numericalAllowance=GUARD)
