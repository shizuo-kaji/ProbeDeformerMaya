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

 
MTypeId probeDeformerNode::id( 0x00000103 );
MString probeDeformerNode::nodeName( "probe" );
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
//    CHECK_MSTATUS_AND_RETURN_IT( status );
    MArrayDataHandle hMatrixArray = data.inputArrayValue(aMatrix);
    MArrayDataHandle hInitMatrixArray = data.inputArrayValue(aInitMatrix);
    int num = hMatrixArray.elementCount();
    if(num != hInitMatrixArray.elementCount() || num==0 || blendMode == 99)
        return MS::kSuccess;
	if( ! rotationCosistency || num != prevNs.size())
	{
		prevThetas.clear();
		prevThetas.resize(num, 0.0);
		prevNs.clear();
		prevNs.resize(num, Vector3d::Zero());
	}
// setting transformation matrix
    std::vector<Matrix4d> initMatrix(num), matrix(num);
    readMatrixArray(hInitMatrixArray, initMatrix);
    readMatrixArray(hMatrixArray, matrix);
    std::vector<Matrix4d> aff(num);
    std::vector<Vector3d> center(num);
    for( int i=0;i<num;i++){
        aff[i]=initMatrix[i].inverse()*matrix[i];
        center[i] << transPart(initMatrix[i]);
    }
    std::vector<Matrix3d> logR(num),R(num),logS(num),logGL(num);
    std::vector<Matrix4d> logSE(num),SE(num),logAff(num);
    std::vector<Vector3d> L(num);
    std::vector<Vector4d> quat(num);
    if(blendMode == 0 || blendMode == 1 || blendMode == 5)  // polarexp or quaternion
    {
        for( int i=0;i<num;i++){
            parametriseGL(aff[i].block(0,0,3,3), logS[i] ,R[i]);
            L[i] = transPart(aff[i]);
            if(blendMode == 0){  // Rotational log
                logR[i]=logSOc(R[i], prevThetas[i], prevNs[i]);
            }else if(blendMode == 1){ // Eucledian log
                SE[i]=affine(R[i], L[i]);
                logSE[i]=logSEc(SE[i], prevThetas[i], prevNs[i]);
            }else if(blendMode == 5){ // quaternion
                Quaternion<double> Q(R[i].transpose());
                quat[i] << Q.x(), Q.y(), Q.z(), Q.w();
            }
        }
    }else if(blendMode == 2){    //logmatrix3
        for( int i=0;i<num;i++){
            logGL[i] = aff[i].block(0,0,3,3).log();
            L[i] = transPart(aff[i]);
        }
    }else if(blendMode == 3){   // logmatrix4
        for( int i=0;i<num;i++){
            logAff[i] = aff[i].log();
        }
    }

    
// transform target vertices
    // get positions
    MPointArray pts;
    itGeo.allPositions(pts);
    int numPts = pts.length();
    if(worldMode){
        for(int j=0; j<numPts; j++ )
            pts[j] *= localToWorldMatrix;
    }
#pragma omp parallel for
    for(int j=0; j<numPts; j++ ){
        // weight computation
        double sidist=0;
        std::vector<double> idist(num);
		Vector3d p;
        p << pts[j].x, pts[j].y, pts[j].z;

        for( int i=0; i<num; i++){
            idist[i] = 1.0 / pow((p-center[i]).squaredNorm(),normExponent/2.0);
            sidist += idist[i];
        }
        std::vector<double> wr(num),ws(num),wl(num);
        if(weightMode == 0){
            for( int i=0; i<num; i++){
                wr[i] = ws[i] = wl[i] = idist[i]/sidist;
            }
        }else{
            float val;
            for( int i=0; i<num; i++){
                rWeightCurveR.getValueAtPosition((1.0/(sqrt(idist[i])*maxDist)), val );
                wr[i] = val;
                rWeightCurveS.getValueAtPosition((1.0/(sqrt(idist[i])*maxDist)), val );
                ws[i] = val;
                rWeightCurveL.getValueAtPosition((1.0/(sqrt(idist[i])*maxDist)), val );
                wl[i] = val;
            }
        }
        // blend matrix
		Matrix4d mat=Matrix4d::Zero();
        if(blendMode==0){
            Matrix3d RR=Matrix3d::Zero();
            Matrix3d SS=Matrix3d::Zero();
            Vector3d l=Vector3d::Zero();
            for( int i=0; i<num; i++){
                RR += wr[i] * logR[i];
                SS += ws[i] * logS[i];
                l += wl[i] * L[i];
            }
            SS = expSym(SS);
            if(frechetSum){
                RR = frechetSO(R, wr);
            }else{
                RR = expSO(RR);
            }
            mat = affine(SS*RR, l);
        }else if(blendMode==1){    // rigid transformation
            Matrix4d EE=Matrix4d::Zero();
            Matrix3d SS=Matrix3d::Zero();
            for( int i=0; i<num; i++){
                EE +=  wr[i] * logSE[i];
                SS +=  ws[i] * logS[i];
            }
            if(frechetSum){
                EE = frechetSE(SE, wr);
            }else{
                EE = expSE(EE);
            }
            mat = affine(expSym(SS),Vector3d::Zero())*EE;
        }else if(blendMode == 2){    //logmatrix3
            Matrix3d G=Matrix3d::Zero();
            Vector3d l=Vector3d::Zero();
            for( int i=0; i<num; i++){
                G +=  wr[i] * logGL[i];
                l += wl[i] * L[i];
            }
            mat = affine(G.exp(), l);
        }else if(blendMode == 3){   // logmatrix4
            Matrix4d A=Matrix4d::Zero();
            for( int i=0; i<num; i++)
                A +=  wr[i] * logAff[i];
            mat = A.exp();
        }else if(blendMode == 5){ // quaternion
            Vector4d q=Vector4d::Zero();
            Matrix3d SS=Matrix3d::Zero();
            Vector3d l=Vector3d::Zero();
            for( int i=0; i<num; i++){
                q += wr[i] * quat[i];
                SS += ws[i] * logS[i];
                l += wl[i] * L[i];
            }
            SS = expSym(SS);
            Quaternion<double> Q(q);
            Matrix3d RR = Q.matrix().transpose();
            mat = affine(SS*RR, l);
        }else if(blendMode==10){
            for( int i=0; i<num; i++){
                mat += wr[i] * aff[i];
            }
        }
        // apply matrix
        MMatrix m;
		for( int i=0; i<4; i++){
			for( int k=0; k<4; k++){
				m(i,k)=mat(i,k);
			}
		}
        pts[j] *= m;
        if(worldMode)
            pts[j] *= localToWorldMatrix.inverse();
    }
    // set positions
    itGeo.setAllPositions(pts);
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

    aBlendMode = eAttr.create( "blendMode", "bm", 0 );
    eAttr.addField( "polarexp", 0 );
    eAttr.addField( "polarexpSE", 1 );
    eAttr.addField( "logmatrix3", 2 );
    eAttr.addField( "logmatrix4", 3 );
    eAttr.addField( "quaternion", 5 );
    eAttr.addField( "linear", 10 );
    eAttr.addField( "off", 99 );
    eAttr.setStorable(true);
    addAttribute( aBlendMode );
    attributeAffects( aBlendMode, outputGeom );

	aRotationConsistency = nAttr.create( "rotationConsistency", "rc", MFnNumericData::kBoolean, 0 );
    nAttr.setStorable(true);
    addAttribute( aRotationConsistency );
    attributeAffects( aRotationConsistency, outputGeom );

	aFrechetSum = nAttr.create( "frechetSum", "fs", MFnNumericData::kBoolean, 0 );
    nAttr.setStorable(true);
    addAttribute( aFrechetSum );
    attributeAffects( aFrechetSum, outputGeom );

	aWorldMode = nAttr.create( "worldMode", "wrldmd", MFnNumericData::kBoolean, 0 );
    nAttr.setStorable(true);
    addAttribute( aWorldMode );
    attributeAffects( aWorldMode, outputGeom );
    
    aWeightMode = eAttr.create( "weightMode", "wtm", 0 );
    eAttr.addField( "normal", 0 );
    eAttr.addField( "curve", 1 );
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

//    aBlendWeight = nAttr.create( "blendWeight", "bw", MFnNumericData::kDouble );
//    nAttr.setKeyable( true );
//    addAttribute( aBlendWeight );
//    attributeAffects( aBlendWeight, outputGeom );
 
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
