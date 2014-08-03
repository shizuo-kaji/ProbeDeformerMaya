/**
 * @file probeDeformer.cpp
 * @brief Probe Deformer plugin for Maya
 * @section LICENSE The MIT License
 * @section  requirements:  Eigen library, Maya
 * @version 0.14
 * @date  3/Nov/2013
 * @author Shizuo KAJI
 */

#include "StdAfx.h"
#include "probeDeformer.h"

using namespace Eigen;
using namespace AffineLib;

// parametrisation mode
#define BM_SRL 0
#define BM_SES 1
#define BM_LOG3 3
#define BM_LOG4 4
#define BM_QSL 5
#define BM_AFF 10
#define BM_OFF -1

// weight mode
#define WM_INV_DISTANCE 0
#define WM_CUTOFF_DISTANCE 1
#define WM_DRAW 2

MTypeId probeDeformerNode::id( 0x00000103 );
MString probeDeformerNode::nodeName( "probeDeformer" );
MObject probeDeformerNode::aMatrix;
MObject probeDeformerNode::aInitMatrix;
MObject probeDeformerNode::aWorldMode;
MObject probeDeformerNode::aBlendMode;
MObject probeDeformerNode::aWeightMode;
MObject probeDeformerNode::aWeightCurveR;
MObject probeDeformerNode::aWeightCurveS;
MObject probeDeformerNode::aWeightCurveL;
MObject probeDeformerNode::aMaxDist;
MObject probeDeformerNode::aRotationConsistency;
MObject probeDeformerNode::aFrechetSum;
MObject probeDeformerNode::aNormExponent;

void* probeDeformerNode::creator() { return new probeDeformerNode; }
 
MStatus probeDeformerNode::deform( MDataBlock& data, MItGeometry& itGeo, const MMatrix &localToWorldMatrix, unsigned int mIndex )
{
	MObject thisNode = thisMObject();
    MStatus status;
    MThreadUtils::syncNumOpenMPThreads();    // for OpenMP
    bool worldMode = data.inputValue( aWorldMode ).asBool();
    short blendMode = data.inputValue( aBlendMode ).asShort();
    short weightMode = data.inputValue( aWeightMode ).asShort();
    double maxDist = data.inputValue( aMaxDist ).asDouble();
    double normExponent = data.inputValue( aNormExponent ).asDouble();
	bool rotationCosistency = data.inputValue( aRotationConsistency ).asBool();
	bool frechetSum = data.inputValue( aFrechetSum ).asBool();
	MRampAttribute rWeightCurveR( thisNode, aWeightCurveR, &status );
	MRampAttribute rWeightCurveS( thisNode, aWeightCurveS, &status );
	MRampAttribute rWeightCurveL( thisNode, aWeightCurveL, &status );
    MArrayDataHandle hMatrixArray = data.inputArrayValue(aMatrix);
    MArrayDataHandle hInitMatrixArray = data.inputArrayValue(aInitMatrix);
    int numPrb = hMatrixArray.elementCount();
    if(numPrb != hInitMatrixArray.elementCount() || numPrb==0 || blendMode == BM_OFF)
        return MS::kSuccess;
	if( ! rotationCosistency || numPrb != prevNs.size()){
		prevThetas.clear();
		prevThetas.resize(numPrb, 0.0);
		prevNs.clear();
		prevNs.resize(numPrb, Vector3d::Zero());
	}
// setting transformation matrix
    std::vector<Matrix4d> initMatrix(numPrb), matrix(numPrb);
    readMatrixArray(hInitMatrixArray, initMatrix);
    readMatrixArray(hMatrixArray, matrix);
    std::vector<Matrix3d> logR(numPrb),R(numPrb),logS(numPrb),S(numPrb),logGL(numPrb);
    std::vector<Matrix4d> logSE(numPrb),SE(numPrb),logAff(numPrb),Aff(numPrb);
    std::vector<Vector3d> L(numPrb);
    std::vector<Vector4d> quat(numPrb);
    std::vector<Vector3d> probeCenter(numPrb);
    for( int i=0;i<numPrb;i++){
        Aff[i]=initMatrix[i].inverse()*matrix[i];
        probeCenter[i] << transPart(initMatrix[i]);
    }
    if(blendMode == BM_SRL || blendMode == BM_SES || blendMode == BM_QSL){
        for(int i=0;i<numPrb;i++){
            parametriseGL(Aff[i].block(0,0,3,3), logS[i] ,R[i]);
            L[i] = transPart(Aff[i]);
            if(blendMode == BM_SRL){
                logR[i]=logSOc(R[i], prevThetas[i], prevNs[i]);
            }else if(blendMode == BM_SES){
                SE[i]=pad(R[i], L[i]);
                logSE[i]=logSEc(SE[i], prevThetas[i], prevNs[i]);
            }else if(blendMode == BM_QSL){
                Quaternion<double> Q(R[i].transpose());
                quat[i] << Q.x(), Q.y(), Q.z(), Q.w();
                S[i]=expSym(logS[i]);
            }
        }
    }else if(blendMode == BM_LOG3){
        for(int i=0;i<numPrb;i++){
            logGL[i] = Aff[i].block(0,0,3,3).log();
            L[i] = transPart(Aff[i]);
        }
    }else if(blendMode == BM_LOG4){
        for(int i=0;i<numPrb;i++){
            logAff[i] = Aff[i].log();
        }
    }

    
// transform target vertices
    // get positions
    MPointArray Mpts;
    itGeo.allPositions(Mpts);
    int numPts = Mpts.length();
    if(worldMode){
        for(int j=0; j<numPts; j++ )
            Mpts[j] *= localToWorldMatrix;
    }
    std::vector<Vector3d> pts(numPts);
    for(int i=0;i<numPts;i++){
        pts[i] << Mpts[i].x, Mpts[i].y, Mpts[i].z;
    }

#pragma omp parallel for
    for(int j=0; j<numPts; j++ ){
        // weight computation
        std::vector<double> dist(numPrb);

        for( int i=0; i<numPrb; i++){
            dist[i] = (pts[j]-probeCenter[i]).norm();
        }
        std::vector<double> wr(numPrb),ws(numPrb),wl(numPrb);
        if(weightMode == WM_INV_DISTANCE){
            double sum=0;
            std::vector<double> idist(numPrb);
            for( int i=0; i<numPrb; i++){
                idist[i] = 1.0/pow(dist[i],normExponent);
                sum += idist[i];
            }
            assert( sum > 0);
            for( int i=0; i<numPrb; i++){
                wr[i] = ws[i] = wl[i] = idist[i]/sum;
            }
        }else if(weightMode == WM_CUTOFF_DISTANCE){
            for( int i=0; i<numPrb; i++){
                wr[i] = ws[i] = wl[i] = (dist[i] > maxDist)
                ? 0 : pow((maxDist-dist[i])/maxDist,normExponent);
            }
        }else if(weightMode == WM_DRAW){
            float val;
            for( int i=0; i<numPrb; i++){
                rWeightCurveR.getValueAtPosition(dist[i]/maxDist, val );
                wr[i] = val;
                rWeightCurveS.getValueAtPosition(dist[i]/maxDist, val );
                ws[i] = val;
                rWeightCurveL.getValueAtPosition(dist[i]/maxDist, val );
                wl[i] = val;
            }
        }
        Matrix4d mat;
        // blend matrix
        if(blendMode == BM_SRL){
            Matrix3d RR,SS=expSym(blendMat(logS, ws));
            Vector3d l=blendMat(L, wl);
            if(frechetSum){
                RR = frechetSO(R, wr);
            }else{
                RR = expSO(blendMat(logR, wr));
            }
            mat = pad(SS*RR, l);
        }else if(blendMode == BM_SES){
            Matrix4d RR;
            Matrix3d SS=expSym(blendMat(logS, ws));
            if(frechetSum){
                RR = frechetSE(SE, wr);
            }else{
                RR = expSE(blendMat(logSE, wr));
            }
            mat = pad(SS,Vector3d::Zero()) * RR;
        }else if(blendMode == BM_LOG3){
            Matrix3d RR=blendMat(logGL, wr).exp();
            Vector3d l=blendMat(L, wl);
            mat = pad(RR, l);
        }else if(blendMode == BM_LOG4){
            mat=blendMat(logAff, wr).exp();
        }else if(blendMode == BM_QSL){
            Vector4d q=blendQuat(quat,wr);
            Vector3d l=blendMat(L, wl);
            Matrix3d SS=blendMatLin(S,ws);
            Quaternion<double> Q(q);
            Matrix3d RR = Q.matrix().transpose();
            mat = pad(SS*RR, l);
        }else if(blendMode == BM_AFF){
            mat = blendMatLin(Aff,wr);
        }
        // apply matrix
        RowVector4d p = pad(pts[j]) * mat;
        Mpts[j].x = p[0];
        Mpts[j].y = p[1];
        Mpts[j].z = p[2];
        if(worldMode)
            Mpts[j] *= localToWorldMatrix.inverse();
    }
    // set positions
    itGeo.setAllPositions(Mpts);
    return MS::kSuccess;
}


// read array of matrix attributes and convert them to Eigen matrices
void probeDeformerNode::readMatrixArray(MArrayDataHandle& handle, std::vector<Matrix4d>& m)
{
    int num=handle.elementCount();
    MMatrix mat;
    for( int i=0;i<num;i++){
        handle.jumpToArrayElement(i);
        mat=handle.inputValue().asMatrix();
        m[i] << mat(0,0), mat(0,1), mat(0,2), mat(0,3),
            mat(1,0), mat(1,1), mat(1,2), mat(1,3),
            mat(2,0), mat(2,1), mat(2,2), mat(2,3),
            mat(3,0), mat(3,1), mat(3,2), mat(3,3);
    }
}

// create nodes
MStatus probeDeformerNode::initialize()
{
    MFnTypedAttribute tAttr;
    MFnNumericAttribute nAttr;
    MFnEnumAttribute eAttr;
    MFnMatrixAttribute mAttr;
   	MRampAttribute rAttr;

    aMatrix = mAttr.create("probeMatrix", "pm");
    mAttr.setStorable(false);
    mAttr.setHidden(true);
    mAttr.setArray(true);
    mAttr.setUsesArrayDataBuilder(true);
    addAttribute(aMatrix);
    attributeAffects( aMatrix, outputGeom );

    aInitMatrix = mAttr.create("initProbeMatrix", "ipm");
    mAttr.setHidden(true);
    mAttr.setArray(true);
    mAttr.setStorable(true);
    mAttr.setUsesArrayDataBuilder(true);
    addAttribute(aInitMatrix);
    attributeAffects( aInitMatrix, outputGeom );

    aBlendMode = eAttr.create( "blendMode", "bm", BM_SRL );
    eAttr.addField( "expSO+expSym", BM_SRL );
    eAttr.addField( "expSE+expSym", BM_SES );
    eAttr.addField( "logmatrix3", BM_LOG3 );
    eAttr.addField( "logmatrix4", BM_LOG4 );
    eAttr.addField( "quat+linear", BM_QSL );
    eAttr.addField( "linear", BM_AFF );
    eAttr.addField( "off", BM_OFF );
    eAttr.setStorable(true);
    addAttribute( aBlendMode );
    attributeAffects( aBlendMode, outputGeom );

	aRotationConsistency = nAttr.create( "rotationConsistency", "rc", MFnNumericData::kBoolean, false );
    nAttr.setStorable(true);
    addAttribute( aRotationConsistency );
    attributeAffects( aRotationConsistency, outputGeom );

	aFrechetSum = nAttr.create( "frechetSum", "fs", MFnNumericData::kBoolean, false );
    nAttr.setStorable(true);
    addAttribute( aFrechetSum );
    attributeAffects( aFrechetSum, outputGeom );

	aWorldMode = nAttr.create( "worldMode", "wrldmd", MFnNumericData::kBoolean, true );
    nAttr.setStorable(true);
    addAttribute( aWorldMode );
    attributeAffects( aWorldMode, outputGeom );
    
    aWeightMode = eAttr.create( "weightMode", "wtm", WM_INV_DISTANCE );
    eAttr.addField( "inverse", WM_INV_DISTANCE );
    eAttr.addField( "cutoff", WM_CUTOFF_DISTANCE );
    eAttr.addField( "draw", WM_DRAW );
    eAttr.setStorable(true);
    addAttribute( aWeightMode );
    attributeAffects( aWeightMode, outputGeom );

	aMaxDist = nAttr.create("maxDistance", "md", MFnNumericData::kDouble, 8.0);
    nAttr.setMin( 0.001 );
    nAttr.setStorable(true);
	addAttribute( aMaxDist );
	attributeAffects( aMaxDist, outputGeom );

	aNormExponent = nAttr.create("normExponent", "ne", MFnNumericData::kDouble, 1.0);
    nAttr.setStorable(true);
	addAttribute( aNormExponent );
	attributeAffects( aNormExponent, outputGeom );
    
	//ramp
    aWeightCurveR = rAttr.createCurveRamp( "weightCurveRotation", "wcr" );
    addAttribute( aWeightCurveR );
	attributeAffects( aWeightCurveR, outputGeom );
    aWeightCurveS = rAttr.createCurveRamp( "weightCurveShear", "wcs" );
    addAttribute( aWeightCurveS );
	attributeAffects( aWeightCurveS, outputGeom );
    aWeightCurveL = rAttr.createCurveRamp( "weightCurveTranslation", "wcl" );
    addAttribute( aWeightCurveL );
	attributeAffects( aWeightCurveL, outputGeom );

    return MS::kSuccess;
}

// create ramp attributes
MStatus probeDeformerNode::accessoryNodeSetup(MDagModifier& cmd)
{
	MStatus stat;
	MObject thisNode = thisMObject();
	MRampAttribute rWeightCurveR( thisNode, probeDeformerNode::aWeightCurveR, &stat );
	MRampAttribute rWeightCurveS( thisNode, probeDeformerNode::aWeightCurveS, &stat );
	MRampAttribute rWeightCurveL( thisNode, probeDeformerNode::aWeightCurveL, &stat );
	
	MFloatArray a1,b1;// position, value
	MIntArray c1;// interpolation

	a1.append(float(0.0));
    a1.append(float(1.0));
    
    b1.append(float(1.0));
    b1.append(float(0.0));
    
	c1.append(MRampAttribute::kSmooth);
	c1.append(MRampAttribute::kSmooth);
    
    rWeightCurveR.addEntries(a1,b1,c1);
    rWeightCurveS.addEntries(a1,b1,c1);
    rWeightCurveL.addEntries(a1,b1,c1);
    return stat;
}



// plugin initializer
MStatus initializePlugin( MObject obj )
{
    MStatus status;
    MFnPlugin plugin( obj, "Shizuo KAJI", "0.1", "Any");
    status = plugin.registerNode( probeDeformerNode::nodeName, probeDeformerNode::id, probeDeformerNode::creator, probeDeformerNode::initialize, MPxNode::kDeformerNode );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    return status;
}
 
MStatus uninitializePlugin( MObject obj )
{
    MStatus   status;
    MFnPlugin plugin( obj );
    status = plugin.deregisterNode( probeDeformerNode::id );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    return status;
}
