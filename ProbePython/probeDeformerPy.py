# -*- coding: utf-8 -*-

#  Python version of "probeDeformer" Maya plugin
#  @author      Shizuo KAJI
#  @date        2013/5/13

# from script editor in Maya
# import maya.cmds as cmds
# cmds.createNode("probeDeformer")

import numpy as np
import scipy.linalg as linalg
import lib.matlib as ml
import lib.blendlib as bl
import lib.mayaattr as attr

import maya.OpenMayaMPx as OpenMayaMPx
import maya.OpenMaya as OpenMaya
import maya.cmds as cmds

# Uncomment for debug
#import lib.debugmaya as debugmaya
#debugmaya.startDebug()

class ProbeDeformerNode(OpenMayaMPx.MPxDeformerNode):
    kPluginNodeId = OpenMaya.MTypeId(0x00000003)
    kPluginNodeTypeName = "probeDeformerPy"
     
    aMatrix = OpenMaya.MObject()
    aInitMatrix = OpenMaya.MObject()
    aWeights = OpenMaya.MObject()
    aTransMode = OpenMaya.MObject()
    aBlendMode  = OpenMaya.MObject()
    aRotationConsistency = OpenMaya.MObject()
    aWeightMode = OpenMaya.MObject()
    aWeightCurveR = OpenMaya.MObject()
    aWeightCurveS = OpenMaya.MObject()
    aWeightCurveL = OpenMaya.MObject()
    aMaxDist = OpenMaya.MObject()
    aNormExponent = OpenMaya.MObject()
    
    prevtheta=[]
    prevn=[]
    
    def __init__(self):
        OpenMayaMPx.MPxDeformerNode.__init__(self)
 
    def deform(self, data, itGeo, localToWorldMatrix, mIndex):
        envelope = OpenMayaMPx.cvar.MPxDeformerNode_envelope
        env = data.inputValue(envelope).asFloat()
        transMode = data.inputValue( ProbeDeformerNode.aTransMode ).asShort()
        blendMode = data.inputValue( ProbeDeformerNode.aBlendMode ).asShort()
        weightMode = data.inputValue( ProbeDeformerNode.aWeightMode ).asShort()
        maxDist = data.inputValue( ProbeDeformerNode.aMaxDist ).asFloat()
        normExponent = data.inputValue( ProbeDeformerNode.aNormExponent ).asFloat()
        rotationConsistency = data.inputValue( ProbeDeformerNode.aRotationConsistency ).asShort()
        hMatrixArray = data.inputArrayValue(ProbeDeformerNode.aMatrix)
        hInitMatrixArray = data.inputArrayValue(ProbeDeformerNode.aInitMatrix)
        num = hMatrixArray.elementCount()
        if num != hInitMatrixArray.elementCount() or num==0 or blendMode == 99:
            return
        # previous rotation: theta and direction
        if rotationConsistency == 0 or num != len(self.prevn):
            self.prevtheta = [0.0 for i in range(num)]
            self.prevn = [np.array([0.0,0.0,0.0]) for i in range(num)]
            
        # setting transformation matrix
        Aff = []
        P=attr.readMatrixArray(data, ProbeDeformerNode.aInitMatrix)
        Q=attr.readMatrixArray(data, ProbeDeformerNode.aMatrix)
        logS = []
        l = []
        center = []
        R = [0  for i in range(num)]
        for i in range(num):    
            Aff.append(np.dot(np.linalg.inv(P[i]),Q[i]))
            U, s, R[i] = ml.polar(Aff[i][:3,:3])
            logS.append(U*np.mat(np.diag(np.log(s)))*U.T)
            center.append(P[i][3,:3])

        
        # translation part
        if transMode==0:
            l=[Aff[i][3,:3] for i in range(num)]
        else:
            l=[Q[i][3,:3]-P[i][3,:3] for i in range(num)]
            
        # blend mode
        if blendMode == 0:  # polarexp
            if rotationConsistency == 0:
                logR = [ml.logRot(R[i]) for i in range(num)]
            else:
                logR = [0  for i in range(num)]
                for i in range(num):
                    logR[i], self.prevtheta[i], self.prevn[i] = ml.logRotc(R[i], self.prevtheta[i], self.prevn[i]) 
        elif blendMode == 5:   # quaternion3
            Rq = [ml.rot2quat(R[i]) for i in range(num)]
        elif blendMode == 3:  #log matrix4
            logA = [linalg.logm(Aff[i]) for i in range(num)]
        elif blendMode == 6:  #dual quaternion
            dq = [ml.mat2dq(ml.affine(R[i],l[i])) for i in range(num)]
        elif blendMode == 1:  #polarexpSE
            if rotationConsistency == 0:
                logRL = [ml.logRL(ml.affine(R[i],l[i])) for i in range(num)]
            else:
                logRL = [0  for i in range(num)]
                for i in range(num):
                    logRL[i], self.prevtheta[i], self.prevn[i] = ml.logRLc(ml.affine(R[i],l[i]), self.prevtheta[i], self.prevn[i]) 
            
        # weight
        uPtr = OpenMaya.MScriptUtil()
        uPtr.createFromDouble(1.0)
        ptr=uPtr.asFloatPtr()
        hWeightCurveR = OpenMaya.MRampAttribute(self.thisMObject(), self.aWeightCurveR)            
        hWeightCurveS = OpenMaya.MRampAttribute(self.thisMObject(), self.aWeightCurveS)            
        hWeightCurveL = OpenMaya.MRampAttribute(self.thisMObject(), self.aWeightCurveL)            
        # compute target mesh
        #pts = OpenMaya.MPointArray()
        #itGeo.allPositions(pts)
        #for j in range(pts.length()):
        while not itGeo.isDone():
            pt = itGeo.position() * localToWorldMatrix
            v = np.array([pt.x, pt.y, pt.z])
            w = self.weightValue(data, mIndex, itGeo.index())*env
            if weightMode == 0:
                idist=[1.0/ (sum((v-o)**2) ** (normExponent/2.0)) for o in center]
                sidist=sum(idist)
                weightR=[w*d/sidist for d in idist]
                weightS=weightR
                weightL=weightR
            else:
                weightR=[]
                weightS=[]
                weightL=[]
                for o in center:
                    hWeightCurveR.getValueAtPosition(float(np.linalg.norm(v-o)/(maxDist+0.01)), ptr)
                    weightR.append(w*OpenMaya.MScriptUtil().getFloat(ptr))
                    hWeightCurveS.getValueAtPosition(float(np.linalg.norm(v-o)/(maxDist+0.01)), ptr)
                    weightS.append(w*OpenMaya.MScriptUtil().getFloat(ptr))
                    hWeightCurveL.getValueAtPosition(float(np.linalg.norm(v-o)/(maxDist+0.01)), ptr)
                    weightL.append(w*OpenMaya.MScriptUtil().getFloat(ptr))
            #
            if blendMode == 0:
                v=bl.blendS(logS, weightS, v)
                v=bl.blendR(logR, weightR, v)
                v=bl.blendL(l, weightL, v)
            elif blendMode == 5:
                v=bl.blendS(logS, weightS, v)
                v=bl.blendQ(Rq, weightR, v)
                v=bl.blendL(l, weightL, v)
            elif blendMode == 3:
                v=bl.blendA(logA, weightR, v)
            elif blendMode == 6:
                v=bl.blendS(logS, weightS, v)
                v=bl.blendDQ(dq, weightR, v)
            elif blendMode == 1:
                v=bl.blendS(logS, weightS, v)
                v=bl.blendRL(logRL, weightR, v)
            else:
                R=sum([weightR[i]*Aff[i] for i in range(num)])
                v=np.dot(np.insert(v,3,1.0),R)
            
            pt.x, pt.y, pt.z=[v[0],v[1],v[2]]
            pt *= localToWorldMatrix.inverse()
            itGeo.setPosition(pt)
            itGeo.next()
#
        #itGeo.setAllPositions(pts)
        return
    
    #create weight curve attribute
    def accessoryNodeSetup(self, cmd):
        a1 = OpenMaya.MFloatArray()
        b1 = OpenMaya.MFloatArray()
        c1 = OpenMaya.MIntArray()
        
        a1.append(float(0.0))
        a1.append(float(0.5))
        a1.append(float(1.0))
        
        b1.append(float(1.0))
        b1.append(float(0.5))
        b1.append(float(0.0))
        
        c1.append(OpenMaya.MRampAttribute.kSpline)
        c1.append(OpenMaya.MRampAttribute.kSpline)
        c1.append(OpenMaya.MRampAttribute.kSpline)
        
        hWeightCurve = OpenMaya.MRampAttribute(self.thisMObject(), self.aWeightCurveR)       
        hWeightCurve.addEntries(a1,b1,c1)
        hWeightCurve = OpenMaya.MRampAttribute(self.thisMObject(), self.aWeightCurveS)       
        hWeightCurve.addEntries(a1,b1,c1)
        hWeightCurve = OpenMaya.MRampAttribute(self.thisMObject(), self.aWeightCurveL)       
        hWeightCurve.addEntries(a1,b1,c1)
    
     
def creator():
    return OpenMayaMPx.asMPxPtr(ProbeDeformerNode())
 
def initialize():
    outputGeom = OpenMayaMPx.cvar.MPxDeformerNode_outputGeom
    # probe matrix
    mAttr = OpenMaya.MFnMatrixAttribute()
    ProbeDeformerNode.aMatrix = attr.createMatrixArray(mAttr,"probeMatrix", "pm")
    ProbeDeformerNode.addAttribute( ProbeDeformerNode.aMatrix )    
    ProbeDeformerNode.attributeAffects(ProbeDeformerNode.aMatrix, outputGeom)  
    # Initial probe matrix
    mAttr = OpenMaya.MFnMatrixAttribute()
    ProbeDeformerNode.aInitMatrix = attr.createMatrixArray(mAttr, "initProbeMatrix", "ipm")
    ProbeDeformerNode.addAttribute( ProbeDeformerNode.aInitMatrix )    
    ProbeDeformerNode.attributeAffects(ProbeDeformerNode.aInitMatrix, outputGeom)  
    # previous rotation
    nAttr = OpenMaya.MFnNumericAttribute()
    ProbeDeformerNode.aPrevNArray = attr.create3FloatArray(nAttr, 'prevn', 'prn')
    ProbeDeformerNode.addAttribute( ProbeDeformerNode.aPrevNArray )
    ProbeDeformerNode.attributeAffects(ProbeDeformerNode.aPrevNArray, outputGeom)  
    #
    nAttr = OpenMaya.MFnNumericAttribute()
    ProbeDeformerNode.aPrevThetaArray = attr.createFloatArray(nAttr, 'prevtheta', 'prt')
    ProbeDeformerNode.addAttribute( ProbeDeformerNode.aPrevThetaArray )
    ProbeDeformerNode.attributeAffects(ProbeDeformerNode.aPrevThetaArray, outputGeom)  
    # translation mode
    eAttr = OpenMaya.MFnEnumAttribute()
    ProbeDeformerNode.aTransMode = eAttr.create( "transMode", "tm", 0 )
    eAttr.addField( "world", 0 )
    eAttr.addField( "local", 1 )
    eAttr.setStorable( True )
    ProbeDeformerNode.addAttribute( ProbeDeformerNode.aTransMode)    
    ProbeDeformerNode.attributeAffects(ProbeDeformerNode.aTransMode, outputGeom)  
    # interpolation mode 
    eAttr = OpenMaya.MFnEnumAttribute()
    ProbeDeformerNode.aBlendMode = eAttr.create( "blendMode", "bm", 0 )
    eAttr.addField( "polarexp", 0 )
    eAttr.addField( "polarexpSE", 1 )
    eAttr.addField( "logmatrix3", 2 )
    eAttr.addField( "logmatrix4", 3 )
    eAttr.addField( "quaternion", 5 )
    eAttr.addField( "dualQuaternion", 6 )
    eAttr.addField( "linear", 10 )
    eAttr.addField( "off", 99 )
    eAttr.setStorable( True )
    ProbeDeformerNode.addAttribute( ProbeDeformerNode.aBlendMode)    
    ProbeDeformerNode.attributeAffects(ProbeDeformerNode.aBlendMode, outputGeom)  
    # rotation angle consistency mode
    nAttr = OpenMaya.MFnNumericAttribute()
    ProbeDeformerNode.aRotationConsistency = nAttr.create( "rotationConsistency", "rc", OpenMaya.MFnNumericData.kBoolean, 0 )
    nAttr.setStorable( True )
    ProbeDeformerNode.addAttribute( ProbeDeformerNode.aRotationConsistency)    
    ProbeDeformerNode.attributeAffects(ProbeDeformerNode.aRotationConsistency, outputGeom)  
    # weighting mode
    eAttr = OpenMaya.MFnEnumAttribute()
    ProbeDeformerNode.aWeightMode = eAttr.create( "weightMode", "wtm", 0 )
    eAttr.addField( "normal", 0 )
    eAttr.addField( "curve", 1 )
    eAttr.setStorable( True )
    ProbeDeformerNode.addAttribute( ProbeDeformerNode.aWeightMode)    
    ProbeDeformerNode.attributeAffects(ProbeDeformerNode.aWeightMode, outputGeom)  
    # weight curve
    rAttr = OpenMaya.MRampAttribute()    
    ProbeDeformerNode.aWeightCurveR = rAttr.createCurveRamp("weightCurveRotation", "wcr")
    ProbeDeformerNode.addAttribute( ProbeDeformerNode.aWeightCurveR)    
    ProbeDeformerNode.attributeAffects(ProbeDeformerNode.aWeightCurveR, outputGeom)  
    ProbeDeformerNode.aWeightCurveS = rAttr.createCurveRamp("weightCurveShear", "wcs")
    ProbeDeformerNode.addAttribute( ProbeDeformerNode.aWeightCurveS)    
    ProbeDeformerNode.attributeAffects(ProbeDeformerNode.aWeightCurveS, outputGeom)  
    ProbeDeformerNode.aWeightCurveL = rAttr.createCurveRamp("weightCurveTranslation", "wcl")
    ProbeDeformerNode.addAttribute( ProbeDeformerNode.aWeightCurveL)    
    ProbeDeformerNode.attributeAffects(ProbeDeformerNode.aWeightCurveL, outputGeom)  
    # max distance
    nAttr = OpenMaya.MFnNumericAttribute()
    ProbeDeformerNode.aMaxDist = nAttr.create( "maxDistance", "md", OpenMaya.MFnNumericData.kFloat, 10.0)
    ProbeDeformerNode.addAttribute( ProbeDeformerNode.aMaxDist)    
    ProbeDeformerNode.attributeAffects(ProbeDeformerNode.aMaxDist, outputGeom)  
    # norm exponent
    nAttr = OpenMaya.MFnNumericAttribute()
    ProbeDeformerNode.aNormExponent = nAttr.create( "normExponent", "ne", OpenMaya.MFnNumericData.kFloat, 1.0)
    ProbeDeformerNode.addAttribute( ProbeDeformerNode.aNormExponent)
    ProbeDeformerNode.attributeAffects(ProbeDeformerNode.aNormExponent, outputGeom)

def initializePlugin(obj):
    plugin = OpenMayaMPx.MFnPlugin(obj, 'S.Kaji', '1.0', 'Any')
    try:
        plugin.registerNode(ProbeDeformerNode.kPluginNodeTypeName, ProbeDeformerNode.kPluginNodeId, creator, initialize, OpenMayaMPx.MPxNode.kDeformerNode)
    except:
        raise RuntimeError, 'Failed to register node'
 
def uninitializePlugin(obj):
    plugin = OpenMayaMPx.MFnPlugin(obj)
    try:
        plugin.deregisterNode(ProbeDeformerNode.kPluginNodeId)
    except:
        raise RuntimeError, 'Failed to deregister node'
