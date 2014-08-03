# -*- coding: utf-8 -*-
## @package blendlib
#
#　　Affine 行列をブレンドして点に作用させるライブラリ
#  @author      Shizuo KAJI
#  @date        2013/5/21

# Numpy moduleのインポート
import numpy as np
import scipy.linalg as linalg
import lib.matlib as ml

import sys
E=sys.float_info.epsilon*100

## パラメータ
class BlendParameters:
    def __init__(self, p=[], triangles=[]):
        self.p = p   #頂点座標
        self.name = ''   #(メッシュの)名前
        self.weights = []  #各単体ごとのウェイトのリスト
        self.triangles = triangles #表面の三角形のリスト
        self.tet = []   # 単体のリスト
        self.P = [] # 単体ごとの頂点行列
        self.PI = [] # 上の逆行列(initial mesh 用)
        self.Aff = []  #4x4アフィン変換のリスト
        self.R = [] # A(v)=SR(v)+l
        self.S = [] 
        self.l = []
        self.logAff = [] #それぞれのlog
        self.logR = [] 
        self.logS = []
        self.RL = [] # 回転+平行移動=剛体変形部分
        self.quat = [] # 回転に対応する四元数
        self.dq = [] # 剛体変形に対応する双対四元数
    def fTet(self):
        return faceTetMatrix(self.p, self.triangles)
    def setPI(self):
        self.PI=[np.linalg.inv(self.P[i]) for i in range(len(self.P))]
    def setAff(self, bpInitData):
        self.Aff = [np.dot(bpInitData.PI[i],self.P[i]) for i in range(len(self.P))]
    def setL(self):
        self.l=[self.Aff[i][3,:3] for i in range(len(self.Aff))]
    def setPolar(self):
        self.logS=[]
        self.logR=[]
        self.setL()
        for i in range(len(self.Aff)):
                U, s, R = ml.polar(Aff[i][:3,:3])
                self.logS.append(np.dot(np.dot(U,np.diag(np.log(s))),U.T))
                self.logR.append(ml.logRot(R))
        
        
# ブレンドアルゴリズム
def polarexp(logR, logS, l, weight, p ):
    """
        極分解データ A=SR から、R, S 成分ともに exp でブレンド。
        平行移動成分は独立に線型ブレンド。
    @param logR[]: 回転成分の log
    @param logS[]: shear 成分の log
    @param l[]: 平行移動成分 
    @param weight[]: ウェイト 
    @param p: 作用させる点の座標
    @return: 点の座標
    """
    num=len(logR)
    lR=[weight[j]*logR[j] for j in range(num)]
    lS=[weight[j]*logS[j] for j in range(num)]
    ll=sum([weight[j]*l[j] for j in range(num)])
    pos = np.dot(np.dot(p,ml.expSymm(sum(lS))),ml.expRot(sum(lR)))+ll
    return np.squeeze(np.asarray(pos))

def blendR(logR, weight, p):
    num=len(logR)
    lR=[weight[j]*logR[j] for j in range(num)]
    pos = np.dot(p, ml.expRot(sum(lR)))
    return np.squeeze(np.asarray(pos))[:3]

def blendS(logS, weight, p):
    num=len(logS)
    lS=[weight[j]*logS[j] for j in range(num)]
    pos = np.dot(p,ml.expSymm(sum(lS)))
    return np.squeeze(np.asarray(pos))

def blendL(l, weight, p):
    num=len(l)
    ll=sum([weight[j]*l[j] for j in range(num)])
    return np.squeeze(np.asarray(p+ll))

def blendA(logA, weight, p):
    num=len(logA)
    lA=[weight[j]*logA[j]  for j in range(num)]
    pos=np.dot(np.insert(p,3,1.0), linalg.expm(sum(lA)))
    return np.squeeze(np.asarray(pos))[:3]

def blendRL(logRL, weight, p): 
    num=len(logRL)
    lRL=[weight[i]*logRL[i] for i in range(num)]
    pos=np.dot(np.insert(p,3,1.0), ml.expRL(sum(lRL)))
    return np.squeeze(np.asarray(pos))[:3]

def blendQ(quat, weight, p):
    num=len(quat)
    lq=0
    for j in range(num):
        lq+=weight[j]*quat[j]
    uq=lq.unit()
    v=ml.Quaternion((1,p[0],p[1],p[2]))
    pos=uq*v*uq.conj()
    return np.array(pos.v[1:4]) 

def blendDQ(dq, weight, p):
    num=len(dq)
    ldq=0
    for j in range(num):
        ldq+=weight[j]*dq[j]
    v=ml.DualQuaternion(ml.Quaternion((1,0,0,0)),ml.Quaternion((0,p[0],p[1],p[2])))
    udq=ldq.unit()
    pos=udq*v*udq.dconj()
    return np.array(pos.d.v[1:4]) 


def polarexp4(logRL, logS, weight, p):
    lRL=[weight[j]*logRL[j] for j in range(len(logS))]
    lS=[weight[j]*logS[j]  for j in range(len(logS))]
    S=np.identity(4)
    S[:3,:3]=ml.expSymm(sum(lS))
    pos=np.dot(np.dot(np.insert(p,3,1),S),ml.expRL(sum(lRL)))
    return np.squeeze(np.asarray(pos))[:3]

def logmatrix(logA, weight, p ):
    """
        log をとってブレンド(Alexa)
    @param logA[]: 4x4アフィン行列のlog
    @param weight[]: ウェイト 
    @param p: 作用させる点の座標
    @return: 点の座標
    """
    lA=[weight[j]*logA[j]  for j in range(len(logA))]
    A=linalg.expm(sum(lA))
    pos = np.dot(np.insert(p,3,1.0),A)
    return np.squeeze(np.asarray(pos))[:3]

def linear(Aff, weight, p):
    lA=[weight[j]*Aff[j] for j in range(len(Aff))]
    return np.squeeze(np.asarray(np.dot(np.insert(p,3,1),sum(lA))))[:3]

# 単体関係
def faceTetMatrix(pts, triangles):
    """
        メッシュの表面三角形分割から、法線方向を足した四面体の行列を生成
    @param pts[]: 頂点座標
    @param triangles[]: 三角形分割
    @return: 四面体の頂点を縦ベクトルに持つ行列のリスト
    """
    P=[]
    for tr in triangles:
        n=np.cross(pts[tr[1]]-pts[tr[0]],pts[tr[2]]-pts[tr[0]])
        n=n/np.linalg.norm(n)+pts[tr[0]]
        P.append(np.mat(np.column_stack((
                np.array((pts[tr[0]],pts[tr[1]],pts[tr[2]],n)),
                         (1,1,1,1)))))
    return P

def vertexTet(adjList):
    tet =[]
    for i in range(len(adjList)):
        l=len(adjList[i])
        if l==3:
            tet.append([i,adjList[i][0], adjList[i][1], adjList[i][2]])
        else:
            for j in range(l):
                tet.append([i,adjList[i][j], adjList[i][(j+1) % l], adjList[i][(j+2) % l]])
    return tet

if __name__ == "__main__":
    import time

# ARAP アルゴリズム
def arapHI(PI,tet,n=0,alpha=0.001):
    """
    初期単体行列から、ARAP の H の逆行列を求める
    @param PI: 行列
    @param tet: 単体
    @param n: 点の数      
    @param alpha: 平行移動成分の重み 
    @return: Hの逆行列
    """
    if n==0:
        n=tet[-1][-1]+1
    Hlist = [  np.dot( np.dot(PI[i].T, np.diag([1.0,1.0,1.0,alpha**2])),PI[i]) for i in range(len(PI))]
    H = np.zeros((n,n))
    for i in range(len(PI)):
        for j in range(4):
            for k in range(4):
                H[tet[i][j],tet[i][k]] += Hlist[i][j,k]
    return np.linalg.inv(H)
#
def arapG(At, PI, tet, n=0,alpha=0.001):
    """
    理想変形行列から ARAP の G を求める
    @param At: 理想変形行列
    @param PI: 行列
    @param tet: 単体
    @param n: 点の数      
    @param alpha: 平行移動成分の重み 
    @return: Hの逆行列
    """
    if n==0:
        n=tet[-1][-1]+1
    Glist = [ np.dot(np.dot(At[i].T,np.diag([1.0,1.0,1.0,alpha**2])),PI[i])  for i in range(len(tet))]
    G = np.zeros((3,n))
#   G[:3,0]=(1-t)*p[0]+t*q[0]
    for i in range(len(tet)):
        for j in range(4):
            for k in range(3):
                G[k,tet[i][j]] += Glist[i][k,j]
    return G

def arapAt(logR, logS, l, weight):
    """
        ARAP の理想変形行列を、polarexp で求める。
    @param logR[][]: 回転成分の log
    @param logS[][]: shear 成分の log
    @param l[][]: 平行移動成分 
    @param weights[]: ウェイト 
    @return: 理想変形行列
    """
    At=[]
    for i in range(len(logR[0])):
        lR=[weight[j]*logR[j][i] for j in range(len(logR))]
        lS=[weight[j]*logS[j][i] for j in range(len(logR))]
        lA=np.dot(ml.expSymm(sum(lS)),ml.expRot(sum(lR)))
        ll=sum([weight[j]*l[j][i] for j in range(len(logR))])
        At.append(ml.affine(lA,ll))
    return At

def arap(idealMatrix, HI):
    """
    @param pts[]: 頂点座標
    @param triangles[]: 三角形分割
    @return: 四面体の頂点を縦ベクトルに持つ行列のリスト
    """


