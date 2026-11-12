from pathlib import Path
import ctypes,json,csv,time,hashlib
import numpy as np
from scipy.io import loadmat
from scipy.optimize import differential_evolution
HERE=Path(__file__).resolve().parent;ROOT=HERE
lib=ctypes.CDLL(str(HERE/'engine/libbilinear_distance.dylib'));ptr=ctypes.POINTER(ctypes.c_double)
lib.bd_init.argtypes=[ptr,ctypes.c_int,ptr,ctypes.c_int,ptr];lib.bd_distance.argtypes=[ptr,ctypes.c_int,ctypes.c_double,ctypes.c_int,ptr]
def init(x,y,H):
 x=np.ascontiguousarray(x,float);y=np.ascontiguousarray(y,float);H=np.asfortranarray(H,float)
 assert lib.bd_init(x.ctypes.data_as(ptr),len(x),y.ctypes.data_as(ptr),len(y),H.ctypes.data_as(ptr))==0

def distance(P,tol=.001):
 P=np.ascontiguousarray(P,float);out=np.zeros(12)
 assert lib.bd_distance(P.ctypes.data_as(ptr),len(P),tol,300000,out.ctypes.data_as(ptr))==0
 assert out[0]>=0 and out[0]<=out[1] and out[11]==1
 assert abs(np.linalg.norm(out[2:5]-out[5:8])+1e-9-out[1])<1e-8
 return out
