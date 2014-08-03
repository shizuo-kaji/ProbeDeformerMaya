# -*- coding: utf-8 -*-
## @package matlib
#
#　　Affine 行列の操作関連ライブラリ
#  @author      Shizuo KAJI
#  @date        2013/5/7

# Numpy moduleのインポート
import numpy as np

import sys
E=sys.float_info.epsilon*100

# Uncomment for debug
#import lib.debugmaya as debugmaya
#debugmaya.startDebug()

def logSymm(matrix):
    """
    3次正値実対称行列の対数を、対角化によって与える
    @param matrix: 3次実対称行列
    @return: 実対称行列の対数
    """
    eig,P = np.linalg.eigh(matrix)
    return np.dot(np.dot(P, np.diag((np.log(eig[0]),np.log(eig[1]),np.log(eig[2])))), np.transpose(P))

def expSymm(matrix):
    """
    3次対称行列の指数を、対角化によって与える
    @param matrix: 3次実対称行列
    @return: 指数
    """
    eig,P = np.linalg.eigh(matrix)
    return np.dot(np.dot(P, np.diag((np.exp(eig[0]),np.exp(eig[1]),np.exp(eig[2])))), np.transpose(P))
          
def logRot(matrix, epsilon=E):
    """
    　　回転行列の対数を、Rodrigues' formula によって与える
    @param matrix: 3次回転行列
    @return: 対数
    """
    tr=(matrix[0,0]+matrix[1,1]+matrix[2,2]-1.0)/2
    if tr>=1.0:
        theta=0.0
    elif tr<=-1.0:
        theta=np.pi
    else:
        theta=np.arccos(tr)
    if theta<epsilon:
        return np.zeros((3,3))
    elif np.pi-theta<epsilon:
        ref=np.zeros((3,3))
        ref[0,1]=np.pi
        return ref
    M=0.5/np.sin(theta)*(matrix-np.transpose(matrix))
    return theta*M

def logRL(matrix, epsilon=E):
    """
    　　回転+平行移動の4次行列の対数を、Rodrigues' formula によって与える
    @param matrix: 4次行列
    @return: 対数
    """
    m=matrix[:3,:3]
    l=matrix[3,:3]
    X=logRot(m)
    norm2 = X[0,1]**2+X[0,2]**2+X[1,2]**2
    if norm2>epsilon:
        norm = np.sqrt(norm2)
        l=np.dot(l, np.identity(3)-X/2+(2*np.sin(norm)-(1+np.cos(norm))*norm)/(2*norm2*np.sin(norm)) *np.dot(X,X))
        #    l=np.dot(l, np.linalg.inv(np.identity(3)+(1-np.cos(norm))/norm2*X +  (norm-np.sin(norm))/(norm*norm2)* np.dot(X,X)))
    return np.column_stack((np.vstack((X,l)),(0,0,0,0)))


def logRotc(matrix, prevtheta=0.0, prevn=np.array([0.0,0.0,0.0]), epsilon=E):
    """
    　　回転行列の対数を、回転角の consistency を考慮しつつ、Rodrigues' formula によって与える。
    @param matrix: 3次回転行列
    @param prevtheta: 近づけたい回転角
    @param prevn: 近づけたい回転軸(第二成分は負)
    @return: 対数, 回転角, 軸(第二成分は負)
    """
    tr=(matrix[0,0]+matrix[1,1]+matrix[2,2]-1.0)/2
    if tr>=1.0:
        theta=0.0
    elif tr<=-1.0:
        theta=np.pi
    else:
        theta=np.arccos(tr)
    if theta<epsilon or np.pi-theta<epsilon:
        M=np.array([[0,prevn[0],prevn[1]],[-prevn[0],0,prevn[2]],[-prevn[1],-prevn[2],0]])
    else: 
        M=0.5/np.sin(theta)*(matrix-np.transpose(matrix))
    if M[0,1]*prevn[0]+M[0,2]*prevn[1]+M[1,2]*prevn[2]<0:
        M = -M
        theta = -theta
    while prevtheta-theta>np.pi:
        theta += 2*np.pi
    while theta-prevtheta>np.pi:
        theta -= 2*np.pi
    return theta*M, theta, np.array((M[0,1],M[0,2],M[1,2]))

def logRLc(matrix, prevtheta=0.0, prevn=np.array([0.0,0.0,0.0]), epsilon=E):
    """
    　　回転+平行移動の4次行列の対数を、consistent に Rodrigues' formula によって与える。
    @param matrix: 4次行列
    @return: 対数
    """
    m=matrix[:3,:3]
    l=matrix[3,:3]
    X, theta, n=logRotc(m, prevtheta, prevn)
    norm2 = X[0,1]**2+X[0,2]**2+X[1,2]**2
    norm = np.sqrt(norm2)
    if abs(np.sin(norm))>epsilon:
        l=np.dot(l, np.identity(3)-X/2+(2*np.sin(norm)-(1+np.cos(norm))*norm)/(2*norm2*np.sin(norm)) *np.dot(X,X))
    return np.column_stack((np.vstack((X,l)),(0,0,0,0))), theta, n

def expRot(matrix, epsilon=E):
    """
    　　交代行列の指数を、Rodrigues' formula によって与える
    @param matrix: 3次交代行列
    @return: 指数
    """
    norm2 = matrix[0,1]**2+matrix[0,2]**2+matrix[1,2]**2
    if norm2<epsilon:
        return np.identity(3)
    norm = np.sqrt(norm2)
    return np.identity(3)+np.sin(norm)/norm*matrix + (1-np.cos(norm))/norm2 * np.dot(matrix,matrix)

def expRL(matrix, epsilon=E):
    """
    　　回転+平行移動の log の指数を、Rodrigues' formula によって与える
    @param matrix: 4次行列
    @return: 指数
    """
    m=matrix[:3,:3]
    l=matrix[3,:3]
    norm2 = m[0,1]**2+m[0,2]**2+m[1,2]**2
    if norm2<epsilon:
        return affine(np.identity(3),l)
    norm = np.sqrt(norm2)
    R=np.identity(3)+np.sin(norm)/norm*m + (1-np.cos(norm))/norm2 * np.dot(m,m)
    d=np.dot(l, np.identity(3)+(1-np.cos(norm))/norm2*m +  (norm-np.sin(norm))/(norm*norm2)* np.dot(m,m))
    return affine(R,d)

#
def affine(A,l):
    """
        回転成分と平行移動成分からAffine変換行列を生成
    @param A: 実３次行列
    @param l: 実３次ベクトル
    @return: ４斉次アフィン行列
    """
    return np.column_stack((np.vstack((A,l)),(0,0,0,1)))


def polar(matrix):
    """
    正則行列の極分解を求める
    @param matrix: 正則行列
    @return: U,s,R  ,where matrix=Udiag(S)U^{-1}R
    """
    eig, U = np.linalg.eigh(np.dot(matrix,matrix.T))
    s = np.sqrt(eig)
    R = np.dot(np.dot(np.dot(U,np.diag(np.reciprocal(s))),U.T),matrix)
    return U,s,R
    

# 四元数関連
class Quaternion(object):
    def __init__(self, v=(1,0,0,0)):
        self.v = np.array(map(float,v))
        self.r = float(v[0])
        self.i = np.array(map(float,v[1:4]))
    def __pos__(p): 
        return p
    def __neg__(p): 
        return Quaternion(-p.v)
    def __add__(p, q): 
        return Quaternion(p.v+q.v)
    def __iadd__(self,p): 
        return self+p
    def __sub__(p, q): 
        return Quaternion(p.v-q.v)
    def __mul__(p, q): 
        v=p.r*q.i + q.r*p.i + np.cross(p.i, q.i)
        return Quaternion([p.r*q.r - np.dot(p.i, q.i),v[0],v[1],v[2]])
    def __rmul__(p, s):
        return Quaternion(s*p.v)
    def __radd__(p, s):
        return Quaternion([p.v[0]+s,p.v[1],p.v[2],p.v[3]])
    def __div__(p,q):
        return(p*q.inv())
    def inv(p): 
        return Quaternion(p.conj().v / p.abs2())
    def abs2(p):
        return sum(p.v**2)
    def __abs__(p):
        return np.sqrt(p.abs2())
    def conj(p):
        return Quaternion([p.r, -p.i[0],-p.i[1],-p.i[2]])
    def unit(p):
        n=abs(p)
        if n>0.001:
            return 1/n*p
        else:
            return Quaternion((1,0,0,0))
        
    def __repr__(p):
        return "[%f; %f, %f, %f]" % (p.r, p.i[0], p.i[1], p.i[2])
    def x(self): return self.v[1]
    def y(self): return self.v[2]
    def z(self): return self.v[3]
    def mat(self): return quat2rot(self.v)

class DualQuaternion(object):
    # Quaternion q,d:  ordinary and dual parts
    def __init__(self, q,d):
        self.q = q
        self.d = d
    def __pos__(p): 
        return p
    def __neg__(p): 
        return DualQuaternion(-p.q, -p.d)
    def __add__(p, q): 
        return DualQuaternion(p.q+q.q, p.d+q.d)
    def __iadd__(self,p): 
        return self+p
    def __radd__(p, s):
        return DualQuaternion(s+p.q,p.d)
    def __sub__(p, q): 
        return DualQuaternion(p.q-q.q, p.d-q.d)
    def __mul__(p, q): 
        d=p.d*q.q+p.q*q.d
        return DualQuaternion(p.q*q.q, d)
    def __rmul__(p, s):
        return DualQuaternion(Quaternion(s*p.q.v),Quaternion(s*p.d.v))
    def __div__(p,q):
        return(p*q.inv())
    def inv(p): 
        n = sum(p.q.v**2)
        if n>0.001:
            return DualQuaternion(Quaternion(p.q.conj().v/n),
                              Quaternion(p.d.conj().v/n-2*np.dot(p.q.v,p.d.v)*p.q.conj().v/(n**2)))
        else:
            return DualQuaternion(Quaternion((1,0,0,0)), Quaternion((0,0,0,0)))
    def __abs__(p):
        n=np.sqrt(sum(p.q.v**2))
        d=Quaternion([np.dot(p.q.v, p.d.v)/n,0,0,0])
        return DualQuaternion(Quaternion([n,0,0,0]),d)
    def unit(p):
        return p*abs(p).inv()
    def conj(p):
        return DualQuaternion(p.q.conj(),p.d.conj())
    def dconj(p):
        return DualQuaternion(p.q.conj(),-p.d.conj())
    def __repr__(p):
        return str(p.q)+"+e"+str(p.d)
    def mat(self): return dq2mat(self)
    
    
def rot2quat(matrix):
    """
        回転行列から４元数を得る
    @param matrix: 回転行列
    @return: 4元数のarray
    """
    w = 2.0*np.sqrt(1.0 + matrix[0,0] + matrix[1,1] + matrix[2,2])
    x = (matrix[1,2] - matrix[2,1]) / w
    y = (matrix[2,0] - matrix[0,2]) / w
    z = (matrix[0,1] - matrix[1,0]) / w
    return Quaternion((w/4.0, x, y, z))
        
def quat2rot(quat):
    """
        単位４元数から回転行列を得る
    @param quat: ４元数
    @return: 回転行列
    """
    x2=2.0*(quat[1]**2)
    y2=2.0*(quat[2]**2)
    z2=2.0*(quat[3]**2)
    xy=2.0*quat[1]*quat[2]
    yz=2.0*quat[2]*quat[3]
    xz=2.0*quat[1]*quat[3]
    wx=2.0*quat[0]*quat[1]
    wy=2.0*quat[0]*quat[2]
    wz=2.0*quat[0]*quat[3]
    return np.array([[1.0-y2-z2, xy+wz, xz-wy],[xy-wz,1.0-x2-z2,yz+wx],[xz+wy,yz-wx,1.0-x2-y2]])

def dq2mat(dq):
    """
        単位双対４元数から回転行列と平行移動成分を得る
    @param dq: 単位双対４元数
    @return: 4x4剛体変換行列
    """
    R = quat2rot(dq.q.v)
    l = (dq.q.inv()*dq.d).v*2
    return affine(R,l[1:4])

def mat2dq(matrix):
    """
        4x4剛体変換行列から単位双対４元数を得る
    @param dq: 4x4剛体変換行列
    @return: 単位双対４元数
    """
    matrix=np.asarray(matrix)
    q = rot2quat(matrix[:3,:3])
    d = Quaternion([0,matrix[3,0]/2,matrix[3,1]/2,matrix[3,2]/2])
    return DualQuaternion(q,d*q)

if __name__ == "__main__":
    import scipy.linalg as linalg
    import time

    I = np.identity(3)
    O = np.zeros((3,3))
    A = np.array([[0,1,2],[-1,0,-3],[-2,3,0.0]])
    B = np.array([[1,2,3],[2,3,4],[3,2,5.0]])
    R = linalg.expm(A)
    RL = affine(I, (1.0,2.0,3.0))
    X=linalg.logm(RL)
#    X=affine(O,(1,2,3))
    repeat = 2000
      
    # exp speed comparison
    time1 = time.clock()
    for i in range(repeat):
#        exp2 = linalg.expm(X)
        exp2 = linalg.logm(RL)
    time2 = time.clock()
    for i in range(repeat):
#        exp3 = expRL(X)
        exp3 = logRL(RL)
    time3 = time.clock()
    print "scipy"
    print time2 - time1
    print(exp2)
    print "matlib"
    print time3 - time2
    print(exp3)

