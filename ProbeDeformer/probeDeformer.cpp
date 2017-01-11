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
using namespace Tetrise;

typedef Triplet<double> T;


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
MObject probeDeformerNode::aAreaWeighted;
MObject probeDeformerNode::aNeighbourWeighting;

void* probeDeformerNode::creator() { return new probeDeformerNode; }
 
MStatus probeDeformerNode::deform( MDataBlock& data, MItGeometry& itGeo, const MMatrix &localToWorldMatrix, unsigned int mIndex ){
    
    MObject thisNode = thisMObject();
    MStatus status;
    MThreadUtils::syncNumOpenMPThreads();    // for OpenMP
    bool worldMode = data.inputValue( aWorldMode ).asBool();
    bool areaWeighted = data.inputValue( aAreaWeighted ).asBool();
    short blendMode = data.inputValue( aBlendMode ).asShort();
    short weightMode = data.inputValue( aWeightMode ).asShort();
    short visualisationMode = data.inputValue( aVisualisationMode ).asShort();
    double effectRadius = data.inputValue( aEffectRadius ).asDouble();
    double normExponent = data.inputValue( aNormExponent ).asDouble();
	B.rotationConsistency = data.inputValue( aRotationConsistency ).asBool();
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
    B.setNum(numPrb);
// setting transformation matrix
    std::vector<Matrix4d> initMatrix(numPrb), matrix(numPrb);
    readMatrixArray(hInitMatrixArray, initMatrix);
    readMatrixArray(hMatrixArray, matrix);
    for( int i=0;i<numPrb;i++){
        B.Aff[i]=initMatrix[i].inverse()*matrix[i];
        B.centre[i] << transPart(initMatrix[i]);
    }
    B.parametrise(blendMode);
    
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
        for(int j=0; j<numPts; j++ ){
            wr[j].resize(numPrb);ws[j].resize(numPrb);wl[j].resize(numPrb);
        }
        D.setNum(numPrb, numPts, 0);
        D.computeDistPts(pts, B.centre);
        D.findClosestPts();
        if(weightMode == WM_INV_DISTANCE){
            for(int j=0; j<numPts; j++ ){
                for( int i=0; i<numPrb; i++){
                    wr[j][i] = ws[j][i] = wl[j][i] = probeRadius[i]/pow(D.distPts[i][j],normExponent);
                }
            }
        }else if(weightMode == WM_CUTOFF_DISTANCE){
            for(int j=0; j<numPts; j++ ){
                for( int i=0; i<numPrb; i++){
                    wr[j][i] = ws[j][i] = wl[j][i] = (D.distPts[i][j] > probeRadius[i])
                    ? 0 : pow((probeRadius[i]-D.distPts[i][j])/probeRadius[i],normExponent);
                }
            }
        }else if(weightMode == WM_DRAW){
            MRampAttribute rWeightCurveR( thisNode, aWeightCurveR, &status );
            MRampAttribute rWeightCurveS( thisNode, aWeightCurveS, &status );
            MRampAttribute rWeightCurveL( thisNode, aWeightCurveL, &status );
            float val;
            for(int j=0; j<numPts; j++ ){
                for( int i=0; i<numPrb; i++){
                    rWeightCurveR.getValueAtPosition(D.distPts[i][j]/probeRadius[i], val );
                    wr[j][i] = val;
                    rWeightCurveS.getValueAtPosition(D.distPts[i][j]/probeRadius[i], val );
                    ws[j][i] = val;
                    rWeightCurveL.getValueAtPosition(D.distPts[i][j]/probeRadius[i], val );
                    wl[j][i] = val;
                }
            }
        }else if(weightMode & WM_HARMONIC){
            makeFaceTet(data, input, inputGeom, mIndex, pts, M.tetList, M.tetMatrix, M.tetWeight);
            M.numTet = (int)M.tetList.size()/4;
            constraint.resize(numPrb);
            // the vertex closest to the probe is given probeWeight
            for(int i=0;i<numPrb;i++){
                constraint[i]=T(i,D.closestPts[i],probeWeight[i]);
            }
            // vertices within effectRadius are given probeWeight
            if( data.inputValue( aNeighbourWeighting ).asBool() ){
                for(int i=0;i<numPrb;i++){
                    for(int j=0;j<numPts;j++){
                        if(D.distPts[i][j]<effectRadius){
                            constraint.push_back(T(i,j,probeWeight[i]));
                        }
                    }
                }
            }
            // set boundary condition for weight computation
            int numConstraint=constraint.size();
            M.constraintWeight.resize(numConstraint);
            M.constraintVal.resize(numConstraint,numPrb);
            M.constraintVal.setZero();
            for(int i=0;i<numConstraint;i++){
                M.constraintVal(i,constraint[i].row())=constraint[i].value();
                M.constraintWeight[i] = std::make_pair(constraint[i].col(), constraint[i].value());
            }
            // clear tetWeight
            if(!areaWeighted){
                M.tetWeight.clear();
                M.tetWeight.resize(M.numTet,1.0);
            }
            // solve the laplace equation
            int isError;
            if( weightMode == WM_HARMONIC_ARAP){
                M.computeTetMatrixInverse();
                M.dim = numPts + M.numTet;
                isError = M.ARAPprecompute();
            }else{
                M.dim = numPts;
                isError = M.cotanPrecompute();
            }
            if(isError>0) return MS::kFailure;
            M.harmonicSolve();
            for(int i=0;i<numPrb;i++){
                for(int j=0;j<numPts;j++){
                    wr[j][i] = ws[j][i] = wl[j][i] = M.Sol.coeff(j,i);
                }
            }
        }
        // normalise weights
        short normaliseWeightMode = data.inputValue( aNormaliseWeight ).asShort();
        for(int j=0;j<numPts;j++){
            D.normaliseWeight(normaliseWeightMode, wr[j]);
            D.normaliseWeight(normaliseWeightMode, ws[j]);
            D.normaliseWeight(normaliseWeightMode, wl[j]);
        }
        
        // END of weight computation
        status = data.setClean(aComputeWeight);
    }
    
    // compute the blended transformations at each mesh point
    if(weightMode == WM_HARMONIC_TRANS){
        // set boundary condition
        int numConstraint=constraint.size();
        M.constraintWeight.resize(numConstraint);
        M.constraintVal.resize(numConstraint,12);
        M.constraintVal.setZero();
        for(int i=0;i<numConstraint;i++){
            M.constraintVal.row(i) << B.logR[i](0,1), B.logR[i](0,2), B.logR[i](1,2),
                B.logS[i](0,0), B.logS[i](0,1), B.logS[i](0,2),
                B.logS[i](1,1), B.logS[i](1,2), B.logS[i](2,2),
                B.L[i](0), B.L[i](1), B.L[i](2);
        }
        M.harmonicSolve();
        //
        Matrix4d mat;
        Matrix3d RR,SS;
        Vector3d l;
        // TODO: implement for different blendMode
        for(int j=0; j<numPts; j++ ){
            RR << 0, M.Sol(j,0), M.Sol(j,1),
                -M.Sol(j,0), 0, M.Sol(j,2),
                -M.Sol(j,1), -M.Sol(j,2), 0;
            SS << M.Sol(j,3), M.Sol(j,4), M.Sol(j,5),
                M.Sol(j,4), M.Sol(j,6), M.Sol(j,7),
                M.Sol(j,5), M.Sol(j,7), M.Sol(j,8);
            l << M.Sol(j,9), M.Sol(j,10), M.Sol(j,11);
            mat = pad(expSym(SS)*expSO(RR), l);
            RowVector4d p = pad(pts[j]) * mat;
            Mpts[j].x = p[0];
            Mpts[j].y = p[1];
            Mpts[j].z = p[2];
            if(worldMode)
                Mpts[j] *= localToWorldMatrix.inverse();
        }
    }else{
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
                Vector3d l=blendMat(B.L, wll);
                SS = expSym(blendMat(B.logS, wss));
                RR = frechetSum ? frechetSO(B.R, wrr) : expSO(blendMat(B.logR, wr[j]));
                mat = pad(SS*RR, l);
            }else if(blendMode == BM_SSE){
                Matrix4d RR;
                Matrix3d SS=expSym(blendMat(B.logS, wss));
                RR = expSE(blendMat(B.logSE, wrr));
                mat = pad(SS,Vector3d::Zero()) * RR;
            }else if(blendMode == BM_LOG3){
                Matrix3d RR=blendMat(B.logGL, wrr).exp();
                Vector3d l=blendMat(B.L, wll);
                mat = pad(RR, l);
            }else if(blendMode == BM_LOG4){
                mat=blendMat(B.logAff, wrr).exp();
            }else if(blendMode == BM_SQL){
                Vector4d q=blendQuat(B.quat,wrr);
                Vector3d l=blendMat(B.L, wll);
                Matrix3d SS=blendMatLin(B.S,wss);
                Quaternion<double> Q(q);
                Matrix3d RR = Q.matrix().transpose();
                mat = pad(SS*RR, l);
            }else if(blendMode == BM_AFF){
                mat = blendMatLin(B.Aff,wrr);
            }
            // apply matrix
            RowVector4d p = pad(pts[j]) * mat;
            Mpts[j].x = p[0];
            Mpts[j].y = p[1];
            Mpts[j].z = p[2];
            if(worldMode)
                Mpts[j] *= localToWorldMatrix.inverse();
        }
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
    
    aNeighbourWeighting = nAttr.create( "neighbourWeighting", "nghbrw", MFnNumericData::kBoolean, false );
    nAttr.setStorable(true);
    addAttribute( aNeighbourWeighting );
    attributeAffects( aNeighbourWeighting, outputGeom );
    attributeAffects( aNeighbourWeighting, aComputeWeight );
    
    aNormaliseWeight = eAttr.create( "normaliseWeight", "nw", NM_LINEAR );
    eAttr.addField( "NONE", NM_NONE );
    eAttr.addField( "Linear",  NM_LINEAR );
    eAttr.addField( "Softmax", NM_SOFTMAX );
    eAttr.setStorable(true);
    addAttribute( aNormaliseWeight );
    attributeAffects( aNormaliseWeight, outputGeom );
    attributeAffects( aNormaliseWeight, aComputeWeight );

    aWeightMode = eAttr.create( "weightMode", "wtm", WM_INV_DISTANCE );
    eAttr.addField( "inverse", WM_INV_DISTANCE );
    eAttr.addField( "cutoff", WM_CUTOFF_DISTANCE );
    eAttr.addField( "draw", WM_DRAW );
//    eAttr.addField( "harmonic-arap", WM_HARMONIC_ARAP);
    eAttr.addField( "harmonic-cotan", WM_HARMONIC_COTAN);
//    eAttr.addField( "harmonic-trans", WM_HARMONIC_TRANS);
    eAttr.setStorable(true);
    addAttribute( aWeightMode );
    attributeAffects( aWeightMode, outputGeom );
    attributeAffects( aWeightMode, aComputeWeight );

    aAreaWeighted = nAttr.create( "areaWeighted", "aw", MFnNumericData::kBoolean, false );
    nAttr.setStorable(true);
    addAttribute( aAreaWeighted );
    attributeAffects( aAreaWeighted, outputGeom );
    attributeAffects( aAreaWeighted, aComputeWeight );

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
