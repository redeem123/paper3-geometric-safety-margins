"""Bilinear patch bounds from its two triangles and a Hausdorff error.

The maximum vertical deviation from the diagonal triangulation is |D|/4,
where D=z00-z10-z01+z11. It scales by1/4 after a two-axis bisection.
"""
import numpy as np

def point_tri(points,a,b,c):
    points=np.atleast_2d(points);ab=b-a;ac=c-a;n=np.cross(ab,ac);nn=np.dot(n,n)
    projection=points-((points-a)@n/nn)[:,None]*n
    v=projection-a;d00=np.dot(ab,ab);d01=np.dot(ab,ac);d11=np.dot(ac,ac)
    den=nn # Gram determinant, computed without cancellation.
    u=(d11*(v@ab)-d01*(v@ac))/den;w=(d00*(v@ac)-d01*(v@ab))/den
    inside=(u>=0)&(w>=0)&(u+w<=1)
    candidates=[projection];dists=[np.where(inside,np.sum((points-projection)**2,axis=1),np.inf)]
    for q,r in [(a,b),(b,c),(c,a)]:
        e=r-q;t=np.clip((points-q)@e/np.dot(e,e),0,1);x=q+t[:,None]*e
        candidates.append(x);dists.append(np.sum((points-x)**2,axis=1))
    ds=np.stack(dists,axis=1);ix=np.argmin(ds,axis=1)
    return ds[np.arange(len(points)),ix],np.stack(candidates,axis=1)[np.arange(len(points)),ix]

def seg_edges(path,q,r):
    """All route segments versus one triangle edge; boundary and interior KKT cases."""
    p=path[:-1];u=np.diff(path,axis=0);v=r-q;w=p-q
    uu=np.sum(u*u,axis=1);vv=np.dot(v,v);n=np.cross(u,v);den=np.sum(n*n,axis=1)
    candidates=[]
    for t in [np.zeros(len(p)),np.ones(len(p))]:
        s=np.clip(np.sum((w+t[:,None]*u)*v,axis=1)/vv,0,1);candidates.append((t,s))
    for s in [np.zeros(len(p)),np.ones(len(p))]:
        t=np.clip(-np.sum((w-s[:,None]*v)*u,axis=1)/uu,0,1);candidates.append((t,s))
    t=np.zeros(len(p));s=np.zeros(len(p))
    np.divide(np.sum(np.cross(v,w)*n,axis=1),den,out=t,where=den>0)
    np.divide(np.sum(np.cross(u,w)*n,axis=1),den,out=s,where=den>0)
    inside=(den>0)&(t>=0)&(t<=1)&(s>=0)&(s<=1)
    candidates.append((t,s));ds=[];qs=[]
    for i,(t,s) in enumerate(candidates):
        surface=q+s[:,None]*v;pp=p+t[:,None]*u;d=np.sum((surface-pp)**2,axis=1)
        if i==4:d=np.where(inside,d,np.inf)
        ds.append(d);qs.append(surface)
    D=np.stack(ds,axis=1);ix=np.argmin(D,axis=1)
    return D[np.arange(len(p)),ix],np.stack(qs,axis=1)[np.arange(len(p)),ix]

def polyline_triangle(path,a,b,c):
    ds,qs=point_tri(path,a,b,c);best=int(np.argmin(ds));d=ds[best];q=qs[best]
    for x,y in [(a,b),(b,c),(c,a)]:
        dd,qq=seg_edges(path,x,y);i=int(np.argmin(dd))
        if dd[i]<d:d=dd[i];q=qq[i]
    # Interior intersection with the triangle plane.
    n=np.cross(b-a,c-a);v=np.diff(path,axis=0);den=v@n;t=np.zeros(len(v))
    np.divide((a-path[:-1])@n,den,out=t,where=den!=0)
    valid=(den!=0)&(t>=0)&(t<=1)
    if np.any(valid):
        p=path[:-1][valid]+t[valid,None]*v[valid];dd,qq=point_tri(p,a,b,c);i=int(np.argmin(dd))
        if dd[i]<d:d=dd[i];q=qq[i]
    return float(np.sqrt(max(0,d))),q

def patch_triangle_bound(path,p):
    x0,x1,y0,y1,z00,z10,z01,z11=p
    a=np.array([x0,y0,z00]);b=np.array([x1,y0,z10]);c=np.array([x1,y1,z11]);e=np.array([x0,y1,z01])
    d,q=polyline_triangle(path,a,b,c);d2,q2=polyline_triangle(path,a,c,e)
    if d2<d:d,q=d2,q2
    error=abs(z00-z10-z01+z11)/4
    return max(0,d-error-1e-8),q
