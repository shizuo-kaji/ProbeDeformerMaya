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
#include <set>
#include "deformerConst.h"

using namespace Eigen;
using namespace AffineLib;
using namespace Tetrise;

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
MObject probeDeformerNode::aEffectRadius;
MObject probeDeformerNode::aRotationConsistency;
MObject probeDeformerNode::aFrechetSum;
MObject probeDeformerNode::aNormExponent;
MObject probeDeformerNode::aProbeWeight;
MObject probeDeformerNode::aVisualisationMode;
MObject probeDeformerNode::aComputeWeight;
MObject probeDeformerNode::aVisualisationMultiplier;
MObject probeDeformerNode::aNormaliseWeight;

void* probeDeformerNode::creator() { return new probeDeformerNode; }
 
MStatus probeDeformerNode::deform( MDataBlock& data, MItGeometry& itGeo, const MMatrix &localToWorldMatrix, unsigned int mIndex ){
    
    MObject thisNode = thisMObject();
    MStatus status;
    MThreadUtils::syncNumOpenMPThreads();    // for OpenMP
    bool worldMode = data.inputValue( aWorldMode ).asBool();
    short blendMode = data.inputValue( aBlendMode ).asShort();
    short weightMode = data.inputValue( aWeightMode ).asShort();
    short visualisationMode = data.inputValue( aVisualisationMode ).asShort();
    double effectRadius = data.inputValue( aEffectRadius ).asDouble();
    double normExponent = data.inputValue( aNormExponent ).asDouble();
	bool rotationCosistency = data.inputValue( aRotationConsistency ).asBool();
	bool frechetSum = data.inputValue( aFrechetSum ).asBool();
    double visualisationMultiplier = data.inputValue(aVisualisationMultiplier).asDouble();
    
    MArrayDataHandle hInitMatrixArray = data.inputArrayValue(aInitMatrix);
    MArrayDataHandle hMatrixArray = data.inputArrayValue(aMatrix);    
    // delete disconnected probe's attr
    if(hMatrixArray.elementCount() > hInitMatrixArray.elementCount() || hMatrixArray.elementCount() == 0 || blendMode == BM_OFF){
        return MS::kSuccess;
    }else if(hMatrixArray.elementCount() < hInitMatrixArray.elementCount()){
        std::set<int> indices;
        for(int i=0;i<hInitMatrixArray.elementCount();i++){
            hInitMatrixArray.jumpToArrayElement(i);
            indices.insert(hInitMatrixArray.elementIndex());
        }
        for(int i=0;i<hMatrixArray.elementCount();i++){
            hMatrixArray.jumpToArrayElement(i);
            indices.erase(hMatrixArray.elementIndex());
        }
        deleteAttr(data, aInitMatrix, indices);
        deleteAttr(data, aProbeWeight, indices);
    }
    // get positions
    MPointArray Mpts;
    itGeo.allPositions(Mpts);
    int new_numPts = Mpts.length();
    if(worldMode){
        for(int j=0; j<new_numPts; j++ )
            Mpts[j] *= localToWorldMatrix;
    }
    std::vector<Vector3d> pts(new_numPts);
    for(int i=0;i<new_numPts;i++){
        pts[i] << Mpts[i].x, Mpts[i].y, Mpts[i].z;
    }

    //
    bool isNumProbeChanged = (numPrb != hMatrixArray.elementCount() || numPts != new_numPts);
    numPrb = hMatrixArray.elementCount();
    numPts = new_numPts;
    //
	if( ! rotationCosistency || numPrb != logSE.size() || numPrb != logR.size()){
		logSE.clear();
		logSE.resize(numPrb, Matrix4d::Zero().eval());
		logR.clear();
		logR.resize(numPrb, Matrix3d::Zero().eval());
    }
// setting transformation matrix
    std::vector<Matrix4d> initMatrix(numPrb), matrix(numPrb);
    readMatrixArray(hInitMatrixArray, initMatrix);
    readMatrixArray(hMatrixArray, matrix);
    std::vector<Matrix3d> R(numPrb),logS(numPrb),S(numPrb),logGL(numPrb);
    std::vector<Matrix4d> SE(numPrb),logAff(numPrb),Aff(numPrb);
    std::vector<Vector3d> L(numPrb);
    std::vector<Vector4d> quat(numPrb);
    std::vector<Vector3d> probeCenter(numPrb);
    for( int i=0;i<numPrb;i++){
        Aff[i]=initMatrix[i].inverse()*matrix[i];
        probeCenter[i] << transPart(initMatrix[i]);
    }
    if(blendMode == BM_SRL || blendMode == BM_SSE || blendMode == BM_SQL){
        for(int i=0;i<numPrb;i++){
            parametriseGL(Aff[i].block(0,0,3,3), logS[i] ,R[i]);
            L[i] = transPart(Aff[i]);
            if(blendMode == BM_SRL){
                logR[i]=logSOc(R[i], logR[i]);
                S[i]=expSym(logS[i]);
            }else if(blendMode == BM_SSE){
                SE[i]=pad(R[i], L[i]);
                logSE[i]=logSEc(SE[i], logSE[i]);
            }else if(blendMode == BM_SQL){
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
    // load per vertex weights
    std::vector<double> ptsWeight(numPts);
    for (int i=0; !itGeo.isDone(); itGeo.next()){
        ptsWeight[i++] = weightValue(data, mIndex, itGeo.index());
    }
    

    // weight computation
    if(!data.isClean(aComputeWeight) || isNumProbeChanged){
        // load probe weights
        std::vector<double> probeWeight(numPrb), probeRadius(numPrb);
        MArrayDataHandle handle = data.inputArrayValue(aProbeWeight);
        if(handle.elementCount() != numPrb){
            MGlobal::displayInfo("# of Probes and probeWeight are different");
            return MS::kFailure;
        }
        for(int i=0;i<numPrb;i++){
            handle.jumpToArrayElement(i);
            probeWeight[i]=handle.inputValue().asDouble();
            probeRadius[i] = probeWeight[i] * effectRadius;
        }
        //
        wr.resize(numPts),ws.resize(numPts),wl.resize(numPts);
        std::vector< std::vector<double> > dist(numPts);
        for(int j=0; j<numPts; j++ ){
            wr[j].resize(numPrb);ws[j].resize(numPrb);wl[j].resize(numPrb);
            dist[j].resize(numPrb);
            for( int i=0; i<numPrb; i++){
                dist[j][i] = (pts[j]-probeCenter[i]).norm();
            }
        }
        if(weightMode == WM_INV_DISTANCE){
            for(int j=0; j<numPts; j++ ){
                for( int i=0; i<numPrb; i++){
                    wr[j][i] = ws[j][i] = wl[j][i] = probeRadius[i]/pow(dist[j][i],normExponent);
                }
            }
        }else if(weightMode == WM_CUTOFF_DISTANCE){
            for(int j=0; j<numPts; j++ ){
                for( int i=0; i<numPrb; i++){
                    wr[j][i] = ws[j][i] = wl[j][i] = (dist[j][i] > probeRadius[i])
                    ? 0 : pow((probeRadius[i]-dist[j][i])/probeRadius[i],normExponent);
                }
            }
        }else if(weightMode == WM_DRAW){
            MRampAttribute rWeightCurveR( thisNode, aWeightCurveR, &status );
            MRampAttribute rWeightCurveS( thisNode, aWeightCurveS, &status );
            MRampAttribute rWeightCurveL( thisNode, aWeightCurveL, &status );
            float val;
            for(int j=0; j<numPts; j++ ){
                for( int i=0; i<numPrb; i++){
                    rWeightCurveR.getValueAtPosition(dist[j][i]/probeRadius[i], val );
                    wr[j][i] = val;
                    rWeightCurveS.getValueAtPosition(dist[j][i]/probeRadius[i], val );
                    ws[j][i] = val;
                    rWeightCurveL.getValueAtPosition(dist[j][i]/probeRadius[i], val );
                    wl[j][i] = val;
                }
            }
        }else if(weightMode == WM_HARMONIC){
            std::vector<int> fList,tList;
            std::vector< std::vector<double> > ptsWeight(numPrb), w_tet(numPrb);
            std::vector<Matrix4d> P;
            int d=makeFaceTet(data, input, inputGeom, mIndex, pts, fList, tList, P);
            std::vector< std::map<int,double> > weightConstraint(numPrb);
            std::vector<double> weightConstraintValue(0);
            for(int i=0;i<numPrb;i++){
                weightConstraint[i].clear();
            }
            if( weightMode == WM_HARMONIC_NEIGHBOUR ){
                for(int i=0;i<numPrb;i++){
                    for(int j=0;j<numPts;j++){
                        if(dist[j][i]<effectRadius){
                            weightConstraint[i][j] = 1;
                            weightConstraintValue.push_back(probeWeight[i]);
                        }
                    }
                }
            }else if( weightMode == WM_HARMONIC){
                std::vector<int> closestPts(numPrb);
                for(int i=0;i<numPrb;i++){
                    weightConstraint[i].clear();
                    closestPts[i] = 0;
                    double min_d = HUGE_VAL;
                    for(int j=0;j<numPts;j++){
                        if( dist[j][i] < min_d){
                            min_d = dist[j][i];
                            closestPts[i] = j;
                        }
                    }
                }
                for(int i=0;i<numPrb;i++){
                    weightConstraint[i][closestPts[i]] = 1;
                    weightConstraintValue.push_back(probeWeight[i]);
                }
            }
            int isError = harmonicWeight(d, P, tList, fList, weightConstraint, weightConstraintValue, ptsWeight);
            if(isError>0) return MS::kFailure;
            for(int i=0;i<numPrb;i++){
                for(int j=0;j<numPts;j++){
                    wr[j][i] = ws[j][i] = wl[j][i] = ptsWeight[i][j];
                }
            }
        }
        // normalise
        bool normaliseWeight = data.inputValue( aNormaliseWeight ).asBool();
        for(int j=0;j<numPts;j++){
            double sum = std::accumulate(wr[j].begin(), wr[j].end(), 0.0);
            if (sum > 1 || normaliseWeight){
                for (int i = 0; i < numPrb; i++){
                    wr[j][i] /= sum;
                }
            }
            sum = std::accumulate(ws[j].begin(), ws[j].end(), 0.0);
            if (sum > 1 || normaliseWeight){
                for (int i = 0; i < numPrb; i++){
                    ws[j][i] /= sum;
                }
            }
            sum = std::accumulate(wl[j].begin(), wl[j].end(), 0.0);
            if (sum > 1 || normaliseWeight){
                for (int i = 0; i < numPrb; i++){
                    wl[j][i] /= sum;
                }
            }
        }
        // END of weight computation
        status = data.setClean(aComputeWeight);
    }
    
#pragma omp parallel for
    for(int j=0; j<numPts; j++ ){
        std::vector<double> wrr(numPrb),wss(numPrb),wll(numPrb);
        for(int i=0;i<numPrb;i++){
            wrr[i]=ptsWeight[j]*wr[j][i];
            wss[i]=ptsWeight[j]*ws[j][i];
            wll[i]=ptsWeight[j]*wl[j][i];
        }
        // blend matrix
        Matrix4d mat;
        if(blendMode == BM_SRL){
            Matrix3d RR,SS;
            Vector3d l=blendMat(L, wll);
            SS = expSym(blendMat(logS, wss));
            RR = frechetSum ? frechetSO(R, wrr) : expSO(blendMat(logR, wr[j]));
            mat = pad(SS*RR, l);
        }else if(blendMode == BM_SSE){
            Matrix4d RR;
            Matrix3d SS=expSym(blendMat(logS, wss));
            RR = expSE(blendMat(logSE, wrr));
            mat = pad(SS,Vector3d::Zero()) * RR;
        }else if(blendMode == BM_LOG3){
            Matrix3d RR=blendMat(logGL, wrr).exp();
            Vector3d l=blendMat(L, wll);
            mat = pad(RR, l);
        }else if(blendMode == BM_LOG4){
            mat=blendMat(logAff, wrr).exp();
        }else if(blendMode == BM_SQL){
            Vector4d q=blendQuat(quat,wrr);
            Vector3d l=blendMat(L, wll);
            Matrix3d SS=blendMatLin(S,wss);
            Quaternion<double> Q(q);
            Matrix3d RR = Q.matrix().transpose();
            mat = pad(SS*RR, l);
        }else if(blendMode == BM_AFF){
            mat = blendMatLin(Aff,wrr);
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
    
    
    // set vertex colour
    if(visualisationMode != VM_OFF){
        std::vector<double> ptsColour(numPts, 0.0);
        if(visualisationMode == VM_STIFFNESS){
            for(int i=0;i<numPts;i++){
                ptsColour[i] = 1.0 - ptsWeight[i];
            }
        }else if(visualisationMode == VM_EFFECT){
            for(int j=0;j<numPts;j++){
//                ptsColour[j] = std::accumulate(wr[j].begin(), wr[j].end(), 0.0);
                ptsColour[j] = visualisationMultiplier * wr[j][numPrb-1];
            }
        }
        visualise(data, outputGeom, ptsColour);
    }

    return MS::kSuccess;
}

// create attr
MStatus probeDeformerNode::initialize()
{
    MFnTypedAttribute tAttr;
    MFnNumericAttribute nAttr;
    MFnEnumAttribute eAttr;
    MFnMatrixAttribute mAttr;
   	MRampAttribute rAttr;

    // this attr will be dirtied when weight recomputation is needed
    aComputeWeight = nAttr.create( "computeWeight", "computeWeight", MFnNumericData::kBoolean, true );
    nAttr.setHidden(true);
    nAttr.setStorable(false);
    nAttr.setKeyable(false);
    addAttribute( aComputeWeight );

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
    eAttr.addField( "expSE+expSym", BM_SSE );
    eAttr.addField( "logmatrix3", BM_LOG3 );
    eAttr.addField( "logmatrix4", BM_LOG4 );
    eAttr.addField( "quat+linear", BM_SQL );
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
    attributeAffects( aWorldMode, aComputeWeight );
    
	aNormaliseWeight = nAttr.create( "normaliseWeight", "nw", MFnNumericData::kBoolean, true );
    nAttr.setStorable(true);
    addAttribute( aNormaliseWeight );
    attributeAffects( aNormaliseWeight, outputGeom );
    attributeAffects( aNormaliseWeight, aComputeWeight );

    aWeightMode = eAttr.create( "weightMode", "wtm", WM_INV_DISTANCE );
    eAttr.addField( "inverse", WM_INV_DISTANCE );
    eAttr.addField( "cutoff", WM_CUTOFF_DISTANCE );
    eAttr.addField( "draw", WM_DRAW );
    eAttr.addField( "harmonic", WM_HARMONIC);
    eAttr.addField( "harmonic-neibour", WM_HARMONIC_NEIGHBOUR);
    eAttr.setStorable(true);
    addAttribute( aWeightMode );
    attributeAffects( aWeightMode, outputGeom );
    attributeAffects( aWeightMode, aComputeWeight );

	aEffectRadius = nAttr.create("effectRadius", "er", MFnNumericData::kDouble, 8.0);
    nAttr.setMin( EPSILON );
    nAttr.setStorable(true);
	addAttribute( aEffectRadius );
	attributeAffects( aEffectRadius, outputGeom );
	attributeAffects( aEffectRadius, aComputeWeight );

	aNormExponent = nAttr.create("normExponent", "ne", MFnNumericData::kDouble, 1.0);
    nAttr.setStorable(true);
	addAttribute( aNormExponent );
	attributeAffects( aNormExponent, outputGeom );
	attributeAffects( aNormExponent, aComputeWeight );
    
    aProbeWeight = nAttr.create("probeWeight", "prw", MFnNumericData::kDouble, 1.0);
    nAttr.setArray(true);
    nAttr.setStorable(true);
    nAttr.setUsesArrayDataBuilder(true);
    addAttribute(aProbeWeight);
	attributeAffects( aProbeWeight, outputGeom );
	attributeAffects( aProbeWeight, aComputeWeight);

    aVisualisationMode = eAttr.create( "visualisationMode", "vm", VM_OFF );
    eAttr.addField( "off", VM_OFF );
    eAttr.addField( "effect", VM_EFFECT );
    eAttr.addField( "stiffness", VM_STIFFNESS );
    eAttr.setStorable(true);
    addAttribute( aVisualisationMode );
    attributeAffects( aVisualisationMode, outputGeom );
    
    aVisualisationMultiplier = nAttr.create("visualisationMultiplier", "vmp", MFnNumericData::kDouble, 1.0);
    nAttr.setStorable(true);
	addAttribute( aVisualisationMultiplier );
	attributeAffects( aVisualisationMultiplier, outputGeom );

	//ramp
    aWeightCurveR = rAttr.createCurveRamp( "weightCurveRotation", "wcr" );
    addAttribute( aWeightCurveR );
	attributeAffects( aWeightCurveR, outputGeom );
	attributeAffects( aWeightCurveR, aComputeWeight );
    aWeightCurveS = rAttr.createCurveRamp( "weightCurveShear", "wcs" );
    addAttribute( aWeightCurveS );
	attributeAffects( aWeightCurveS, outputGeom );
	attributeAffects( aWeightCurveS, aComputeWeight );
    aWeightCurveL = rAttr.createCurveRamp( "weightCurveTranslation", "wcl" );
    addAttribute( aWeightCurveL );
	attributeAffects( aWeightCurveL, outputGeom );
	attributeAffects( aWeightCurveL, aComputeWeight );

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


// this deformer also changes colours
void probeDeformerNode::postConstructor(){
	setDeformationDetails(kDeformsColors);
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
