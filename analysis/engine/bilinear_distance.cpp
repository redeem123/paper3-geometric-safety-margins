// Complete polyline-to-bilinear-graph numerical bounds. No search algorithm here.
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <queue>
#include <stdexcept>
#include <vector>
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif
struct V { double x,y,z; V operator+(V b)const{return{x+b.x,y+b.y,z+b.z};} V operator-(V b)const{return{x-b.x,y-b.y,z-b.z};} V operator*(double t)const{return{x*t,y*t,z*t};} };
static double dot(V a,V b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static V cross(V a,V b){return{a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static double norm(V a){return std::sqrt(std::max(0.,dot(a,a)));}
static double clamp(double t,double a=0,double b=1){return std::max(a,std::min(b,t));}
static double coord(V a,int k){return k==0?a.x:(k==1?a.y:a.z);}
struct Patch{double x0,x1,y0,y1,a,b,c,d;};
struct Box{V lo,hi;};
struct Node{Box box;int l=-1,r=-1;int i0,i1,j0,j1;};
struct Result{double lower,upper;V path,surface;int segment;double fraction;int pops;bool converged;};
static V clip(V p,Box b){return{clamp(p.x,b.lo.x,b.hi.x),clamp(p.y,b.lo.y,b.hi.y),clamp(p.z,b.lo.z,b.hi.z)};}
static double pointSegment(V p,V a,V b,V& q,double& t){V d=b-a;double dd=dot(d,d);t=dd>0?clamp(dot(p-a,d)/dd):0;q=a+d*t;return norm(q-p);}
static double segmentBox(V a,V b,Box box,V& witness){
 V v=b-a;double t[8]={0,1,0,0,0,0,0,0};int k=2;
 for(int c=0;c<3;c++)for(double side:{coord(box.lo,c),coord(box.hi,c)})t[k++]=coord(v,c)!=0?clamp((side-coord(a,c))/coord(v,c)):0;
 std::sort(t,t+8);double best=INFINITY;
 for(int j=0;j<7;j++){
  double mid=(t[j]+t[j+1])/2,aa=0,vv=0;V p=a+v*mid;
  for(int c=0;c<3;c++){double q=coord(p,c),lower=coord(box.lo,c),upper=coord(box.hi,c),u=coord(v,c);if(q<lower||q>upper){aa+=(coord(a,c)-(q<lower?lower:upper))*u;vv+=u*u;}}
  double s=vv>0?clamp(-aa/vv,t[j],t[j+1]):mid;V p2=a+v*s,q=clip(p2,box);double d=norm(p2-q);if(d<best){best=d;witness=q;}
 }return std::max(0.,best-1e-9);
}
static double boxDistance(const std::vector<V>& path,Box b,V& q){double best=INFINITY;for(size_t i=0;i<std::max(size_t(1),path.size()-1);i++){V w;double d=segmentBox(path[i],path[std::min(i+1,path.size()-1)],b,w);if(d<best){best=d;q=w;}}return best;}
static double pointTriangle(V p,V a,V b,V c,V& q){
 V u=b-a,v=c-a,n=cross(u,v);double nn=dot(n,n);V proj=p-n*(dot(p-a,n)/nn),r=proj-a;double uv=dot(u,v),uu=dot(u,u),vv=dot(v,v);
 double s=(vv*dot(r,u)-uv*dot(r,v))/nn,t=(uu*dot(r,v)-uv*dot(r,u))/nn,best=INFINITY;
 if(s>=0&&t>=0&&s+t<=1){best=norm(p-proj);q=proj;}
 for(auto edge:{std::array<V,2>{a,b},std::array<V,2>{b,c},std::array<V,2>{c,a}}){V w;double f,d=pointSegment(p,edge[0],edge[1],w,f);if(d<best){best=d;q=w;}}
 return best;
}
static double segmentEdge(V a,V b,V q,V r,V& witness){
 V u=b-a,v=r-q,w=a-q,n=cross(u,v);double uu=dot(u,u),vv=dot(v,v),den=dot(n,n),best=INFINITY;
 auto test=[&](double t,double s){V x=a+u*t,y=q+v*s;double d=norm(x-y);if(d<best){best=d;witness=y;}};
 for(double t:{0.,1.})test(t,vv>0?clamp(dot(w+u*t,v)/vv):0);
 for(double s:{0.,1.})test(uu>0?clamp(-dot(w-v*s,u)/uu):0,s);
 if(den>0){double t=dot(cross(v,w),n)/den,s=dot(cross(u,w),n)/den;if(t>=0&&t<=1&&s>=0&&s<=1)test(t,s);}
 return best;
}
static double polyTriangle(const std::vector<V>& path,V a,V b,V c,V& witness){
 double best=INFINITY;for(V p:path){V q;double d=pointTriangle(p,a,b,c,q);if(d<best){best=d;witness=q;}}
 V n=cross(b-a,c-a);
 for(size_t i=0;i+1<path.size();i++){
  for(auto edge:{std::array<V,2>{a,b},std::array<V,2>{b,c},std::array<V,2>{c,a}}){V q;double d=segmentEdge(path[i],path[i+1],edge[0],edge[1],q);if(d<best){best=d;witness=q;}}
  V v=path[i+1]-path[i];double den=dot(v,n);if(den!=0){double t=dot(a-path[i],n)/den;if(t>=0&&t<=1){V q;double d=pointTriangle(path[i]+v*t,a,b,c,q);if(d<best){best=d;witness=q;}}}
 }return best;
}
static Box patchBox(Patch p){return{{p.x0,p.y0,std::min({p.a,p.b,p.c,p.d})},{p.x1,p.y1,std::max({p.a,p.b,p.c,p.d})}};}
static double value(Patch p,double x,double y){double u=clamp((x-p.x0)/(p.x1-p.x0)),v=clamp((y-p.y0)/(p.y1-p.y0));return(1-u)*(1-v)*p.a+u*(1-v)*p.b+(1-u)*v*p.c+u*v*p.d;}
static std::array<Patch,4> split(Patch p){double xm=(p.x0+p.x1)/2,ym=(p.y0+p.y1)/2,zb=(p.a+p.b)/2,zt=(p.c+p.d)/2,zl=(p.a+p.c)/2,zr=(p.b+p.d)/2,zc=(p.a+p.b+p.c+p.d)/4;
 return{{{p.x0,xm,p.y0,ym,p.a,zb,zl,zc},{xm,p.x1,p.y0,ym,zb,p.b,zc,zr},{p.x0,xm,ym,p.y1,zl,zc,p.c,zt},{xm,p.x1,ym,p.y1,zc,zr,zt,p.d}}};}
class Surface{
 public:std::vector<double>x,y,h;std::vector<Node>nodes;int nx,ny;
 Surface(const double* X,int nX,const double* Y,int nY,const double* H):x(X,X+nX),y(Y,Y+nY),h(H,H+nX*nY),nx(nX),ny(nY){if(nx<2||ny<2)throw std::runtime_error("Grid dimensions");for(int j=1;j<nx;j++)if(!(x[j]>x[j-1]))throw std::runtime_error("X not increasing");for(int i=1;i<ny;i++)if(!(y[i]>y[i-1]))throw std::runtime_error("Y not increasing");for(double z:h)if(!std::isfinite(z))throw std::runtime_error("Nonfinite grid");build(0,ny-1,0,nx-1);}
 double height(int i,int j)const{return h[i+ny*j];}
 int build(int i0,int i1,int j0,int j1){Node n;n.i0=i0;n.i1=i1;n.j0=j0;n.j1=j1;double lo=INFINITY,hi=-INFINITY;for(int i=i0;i<=i1;i++)for(int j=j0;j<=j1;j++){double z=height(i,j);lo=std::min(lo,z);hi=std::max(hi,z);}n.box={{x[j0],y[i0],lo},{x[j1],y[i1],hi}};int id=nodes.size();nodes.push_back(n);if(i1-i0>1||j1-j0>1){int l,r;if(i1-i0>=j1-j0&&i1-i0>1){int c=(i0+i1)/2;l=build(i0,c,j0,j1);r=build(c,i1,j0,j1);}else{int c=(j0+j1)/2;l=build(i0,i1,j0,c);r=build(i0,i1,c,j1);}nodes[id].l=l;nodes[id].r=r;}return id;}
 Patch patch(int k)const{auto n=nodes[k];return{x[n.j0],x[n.j1],y[n.i0],y[n.i1],height(n.i0,n.j0),height(n.i0,n.j1),height(n.i1,n.j0),height(n.i1,n.j1)};}
 double z(double xx,double yy)const{int j=clamp(std::upper_bound(x.begin(),x.end(),xx)-x.begin()-1,0,nx-2),i=clamp(std::upper_bound(y.begin(),y.end(),yy)-y.begin()-1,0,ny-2);Patch p{x[j],x[j+1],y[i],y[i+1],height(i,j),height(i,j+1),height(i+1,j),height(i+1,j+1)};return value(p,xx,yy);}
 Result distance(const std::vector<V>& path,double tolerance,int maxPops)const{
  if(path.empty())throw std::runtime_error("Empty path");Result result{0,INFINITY,{},{},0,0,0,false};
  auto upper=[&](V surface){for(size_t i=0;i<std::max(size_t(1),path.size()-1);i++){V q;double t,d=pointSegment(surface,path[i],path[std::min(i+1,path.size()-1)],q,t)+1e-9;if(d<result.upper){result.upper=d;result.path=q;result.surface=surface;result.segment=i;result.fraction=t;}}};
  for(V p:path){if(!std::isfinite(p.x)||!std::isfinite(p.y)||!std::isfinite(p.z)||p.x<x.front()||p.x>x.back()||p.y<y.front()||p.y>y.back())throw std::runtime_error("Path outside grid");upper({p.x,p.y,z(p.x,p.y)});}
  struct Task{double lower;int node;Patch p;bool operator<(const Task& b)const{return lower>b.lower;}};
  std::priority_queue<Task>heap;V w;heap.push({boxDistance(path,nodes[0].box,w),0,{}});
  while(!heap.empty()&&result.upper-heap.top().lower>tolerance&&result.pops<maxPops){Task task=heap.top();heap.pop();result.pops++;if(task.lower>=result.upper)continue;
   if(task.node>=0&&nodes[task.node].l>=0){for(int child:{nodes[task.node].l,nodes[task.node].r}){V q;double lb=boxDistance(path,nodes[child].box,q);upper({q.x,q.y,z(q.x,q.y)});if(lb<result.upper)heap.push({lb,child,{}});}continue;}
   Patch p=task.node>=0?patch(task.node):task.p;
   for(Patch pp:split(p)){
    V q,tri;double lower=boxDistance(path,patchBox(pp),q);
    V a{pp.x0,pp.y0,pp.a},b{pp.x1,pp.y0,pp.b},c{pp.x1,pp.y1,pp.d},e{pp.x0,pp.y1,pp.c},q1,q2;
    double d1=polyTriangle(path,a,b,c,q1),d2=polyTriangle(path,a,c,e,q2);tri=d1<=d2?q1:q2;
    lower=std::max(lower,std::max(0.,std::min(d1,d2)-std::abs(pp.a-pp.b-pp.c+pp.d)/4-1e-8));
    upper({tri.x,tri.y,value(pp,tri.x,tri.y)});upper({q.x,q.y,value(pp,q.x,q.y)});
    upper({(pp.x0+pp.x1)/2,(pp.y0+pp.y1)/2,(pp.a+pp.b+pp.c+pp.d)/4});for(V corner:{a,b,c,e})upper(corner);
    if(lower<result.upper)heap.push({lower,-1,pp});
   }
  }
  result.lower=heap.empty()?result.upper:std::min(result.upper,heap.top().lower);result.converged=result.upper-result.lower<=tolerance;return result;
 }
};
static std::unique_ptr<Surface> surface;
static void flatten(const Result&r,double*out){double a[]={r.lower,r.upper,r.path.x,r.path.y,r.path.z,r.surface.x,r.surface.y,r.surface.z,double(r.segment),r.fraction,double(r.pops),double(r.converged)};std::copy(a,a+12,out);}
extern "C" int bd_init(const double*x,int nx,const double*y,int ny,const double*h){try{surface.reset(new Surface(x,nx,y,ny,h));return 0;}catch(...){return 1;}}
extern "C" int bd_distance(const double*p,int n,double tol,int maxPops,double*out){try{if(!surface||n<1||tol<=0||maxPops<1)return 2;std::vector<V>path;for(int i=0;i<n;i++)path.push_back({p[3*i],p[3*i+1],p[3*i+2]});flatten(surface->distance(path,tol,maxPops),out);return 0;}catch(...){return 1;}}
#ifdef MATLAB_MEX_FILE
void mexFunction(int nlhs,mxArray*plhs[],int nrhs,const mxArray*prhs[]){try{
 if(nrhs<1||!mxIsChar(prhs[0])||nlhs>1)mexErrMsgIdAndTxt("BD:arguments","Command required; at most one output.");char*cmd=mxArrayToString(prhs[0]);std::string command(cmd);mxFree(cmd);
 if(command=="init"){if(nrhs!=4||!mxIsDouble(prhs[1])||!mxIsDouble(prhs[2])||!mxIsDouble(prhs[3])||mxGetNumberOfElements(prhs[1])!=mxGetN(prhs[3])||mxGetNumberOfElements(prhs[2])!=mxGetM(prhs[3]))mexErrMsgIdAndTxt("BD:grid","Use init,x,y,H.");if(bd_init(mxGetPr(prhs[1]),mxGetNumberOfElements(prhs[1]),mxGetPr(prhs[2]),mxGetNumberOfElements(prhs[2]),mxGetPr(prhs[3])))throw std::runtime_error("Invalid grid");if(nlhs)plhs[0]=mxCreateDoubleScalar(surface->nodes.size());return;}
 if(command!="distance"||nrhs!=4||!mxIsDouble(prhs[1])||mxGetN(prhs[1])!=3)throw std::runtime_error("Use distance,N-by-3 path,tolerance,maxPops");int n=mxGetM(prhs[1]);const double*P=mxGetPr(prhs[1]);std::vector<double>row(3*n);for(int i=0;i<n;i++)for(int j=0;j<3;j++)row[3*i+j]=P[i+n*j];double out[12];if(bd_distance(row.data(),n,mxGetScalar(prhs[2]),mxGetScalar(prhs[3]),out))throw std::runtime_error("Distance failed");if(nlhs){plhs[0]=mxCreateDoubleMatrix(1,12,mxREAL);std::copy(out,out+12,mxGetPr(plhs[0]));}
 }catch(const std::exception&e){mexErrMsgIdAndTxt("BD:failure","%s",e.what());}}
#endif
