/**
 * @file probeDeformerARAP.cpp
 * @brief Probe Deformer ARAP plugin for Maya
 * @section LICENSE The MIT License
 * @section  requirements:  Eigen library, Maya
 * @version 0.15
 * @date  3/Jan/2014
 * @author Shizuo KAJI
 *
 * @section Shape is deformed according to the blended affine transformation and ARAP energy.
 * The closest point on the mesh is constraint by the probe.
 */


#include "StdAfx.h"
#include "probeDeformerARAP.h"

using namespace Eigen;
using namespace AffineLib;

MTypeId probeDeformerARAPNode::id( 0x00000104 );
MString probeDeformerARAPNode::nodeName( "probeDeformerARAP" );
MObject probeDeformerARAPNode::aMatrix;
MObject probeDeformerARAPNode::aInitMatrix;
MObject probeDeformerARAPNode::aWorldMode;
MObject probeDeformerARAPNode::aBlendMode;
MObject probeDeformerARAPNode::aWeightMode;
MObject probeDeformerARAPNode::aWeightCurveR;
MObject probeDeformerARAPNode::aWeightCurveS;
MObject probeDeformerARAPNode::aWeightCurveL;
MObject probeDeformerARAPNode::aMaxDist;
MObject probeDeformerARAPNode::aTransWeight;
MObject probeDeformerARAPNode::aConstraintWeight;
MObject probeDeformerARAPNode::aRotationConsistency;
MObject probeDeformerARAPNode::aFrechetSum;
MObject probeDeformerARAPNode::aNormExponent;

void* probeDeformerARAPNode::creator() { return new probeDeformerARAPNode; }
 
MStatus probeDeformerARAPNode::deform( MDataBlock& data, MItGeometry& itGeo, const MMatrix &localToWorldMatrix, unsigned int mIndex )
{
	MObject thisNode = thisMObject();
    MStatus status;
    MThreadUtils::syncNumOpenMPThreads();    // for OpenMP
    
    bool new_worldMode = data.inputValue( aWorldMode ).asBool();
    short blendMode = data.inputValue( aBlendMode ).asShort();
    double new_transWeight = data.inputValue( aTransWeight ).asDouble();
    double new_constraintWeight = data.inputValue( aConstraintWeight ).asDouble();
    double new_normExponent = data.inputValue( aNormExponent ).asDouble();
    MArrayDataHandle hMatrixArray = data.inputArrayValue(aMatrix);
    MArrayDataHandle hInitMatrixArray = data.inputArrayValue(aInitMatrix);
    int new_numPrb = hMatrixArray.elementCount();
    if(new_numPrb != hInitMatrixArray.elementCount())
        return MS::kSuccess;
    // read matrices
    std::vector<Matrix4d> initMatrix(new_numPrb), matrix(new_numPrb);
    readMatrixArray(hInitMatrixArray, initMatrix);
    readMatrixArray(hMatrixArray, matrix);
    // new probe connection
    if(numPts == 0 || blendMode == 99 || normExponent != new_normExponent || transWeight != new_transWeight
       || new_worldMode != worldMode || constraintWeight != new_constraintWeight || numPrb != new_numPrb)
    {
        numPrb = new_numPrb;
        worldMode = new_worldMode;
        transWeight = new_transWeight;
        constraintWeight = new_constraintWeight;
        normExponent = new_normExponent;
        // probe position
        probeCenter.resize(numPrb);
        for(int i=0;i<numPrb;i++)
            probeCenter[i] << initMatrix[i](3,0),initMatrix[i](3,1),initMatrix[i](3,2);
        // read mesh data
        MArrayDataHandle hInput = data.outputArrayValue( input, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = hInput.jumpToElement( mIndex );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        MObject oInputGeom = hInput.outputValue().child( inputGeom ).asMesh();
        MFnMesh inputMesh(oInputGeom);
        inputMesh.getPoints( pts );
		numPts=pts.length();
        if(worldMode){
            for(int j=0; j<numPts; j++ )
                pts[j] *= localToWorldMatrix;
        }
        // prepare list of facial tetrahedra
        MIntArray count;
        inputMesh.getTriangles( count, triangles );
		numTet=triangles.length()/3;
		std::vector<Matrix4d> P(numTet);
        PI.resize(numTet);
        tetCenter.resize(numTet);
        tetMatrixC(pts, triangles, P, tetCenter);
        // compute distance between probe and vertex
        sidist.resize(numTet);
        idist.resize(numTet);
        for(int j=0;j<numTet;j++){
            sidist[j]=0;
            idist[j].resize(numPrb);
            for(int i=0; i<numPrb; i++){
                idist[j][i] = 1.0 / pow((tetCenter[j]-probeCenter[i]).squaredNorm(),normExponent/2.0);
                sidist[j] += idist[j][i];
            }
        }
        // find constraint points
        if(numTet<numPrb) return(MS::kSuccess);
        constraintTet.resize(numPrb);
        constraintVector.resize(numPrb);
        for(int i=0;i<numPrb;i++){
            constraintTet[i] = 0;
            for(int j=1;j<numTet;j++){
                if(idist[j][i] > idist[constraintTet[i]][i]){
                    constraintTet[i] = j;
                }
            }
            constraintVector[i] << tetCenter[constraintTet[i]](0), tetCenter[constraintTet[i]](1), tetCenter[constraintTet[i]](2), 1.0;
        }
        // precompute arap matrix
		for(int i=0;i<numTet;i++){
            assert(P[i].determinant()!=0);
			PI[i] = P[i].inverse();
        }
        arapHI(PI, triangles);
    }
    // read attributes
    short weightMode = data.inputValue( aWeightMode ).asShort();
    double maxDist = data.inputValue( aMaxDist ).asDouble();
	bool rotationCosistency = data.inputValue( aRotationConsistency ).asBool();
	bool frechetSum = data.inputValue( aFrechetSum ).asBool();
	MRampAttribute rWeightCurveR( thisNode, aWeightCurveR, &status );
	MRampAttribute rWeightCurveS( thisNode, aWeightCurveS, &status );
	MRampAttribute rWeightCurveL( thisNode, aWeightCurveL, &status );
    
    
    
// setting transformation matrix
	if( ! rotationCosistency || numPrb != prevNs.size())
	{
		prevThetas.clear();
		prevThetas.resize(numPrb, 0.0);
		prevNs.clear();
		prevNs.resize(numPrb, Vector3d::Zero());
	}
    std::vector<Matrix4d> aff(numPrb);
    std::vector<Vector3d> center(numPrb);
    for(int i=0;i<numPrb;i++)
        aff[i]=initMatrix[i].inverse()*matrix[i];
    std::vector<Matrix3d> logR(numPrb);
    std::vector<Matrix3d> R(numPrb);
    std::vector<Matrix4d> logSE(numPrb);
    std::vector<Matrix4d> SE(numPrb);
    std::vector<Matrix3d> logS(numPrb);
    std::vector<Vector3d> L(numPrb);
    std::vector<Matrix3d> logGL(numPrb);
    std::vector<Matrix4d> logAff(numPrb);
    std::vector<Vector4d> quat(numPrb);
    if(blendMode == 0 || blendMode == 1 || blendMode == 5)  // polarexp or quaternion
    {
        for(int i=0;i<numPrb;i++)        {
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
        for(int i=0;i<numPrb;i++){
            logGL[i] = aff[i].block(0,0,3,3).log();
            L[i] = transPart(aff[i]);
        }
    }else if(blendMode == 3){   // logmatrix4
        for(int i=0;i<numPrb;i++){
            logAff[i] = aff[i].log();
        }
    }


// prepare transform matrix for each simplex
    std::vector<Matrix4d> At(numTet);
#pragma omp parallel for
    for(int j=0; j<numTet; j++ )
    {
        // weight computation
        std::vector<double> wr(numPrb),ws(numPrb),wl(numPrb);
        if(weightMode == 0){
            for(int i=0; i<numPrb; i++){
                wr[i] = ws[i] = wl[i] = idist[j][i]/sidist[j];
            }
        }else{
            float val;
            for(int i=0; i<numPrb; i++){
                rWeightCurveR.getValueAtPosition((1.0/(sqrt(idist[j][i])*maxDist)), val );
                wr[i] = val;
                rWeightCurveS.getValueAtPosition((1.0/(sqrt(idist[j][i])*maxDist)), val );
                ws[i] = val;
                rWeightCurveL.getValueAtPosition((1.0/(sqrt(idist[j][i])*maxDist)), val );
                wl[i] = val;
            }
        }
        // blend matrix
        if(blendMode==0)
        {
            Matrix3d RR=Matrix3d::Zero();
            Matrix3d SS=Matrix3d::Zero();
            Vector3d l=Vector3d::Zero();
            for(int i=0; i<numPrb; i++){
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
            At[j] = affine(SS*RR, l);
        }else if(blendMode==1){    // rigid transformation
            Matrix4d EE=Matrix4d::Zero();
            Matrix3d SS=Matrix3d::Zero();
            for(int i=0; i<numPrb; i++){
                EE +=  wr[i] * logSE[i];
                SS +=  ws[i] * logS[i];
            }
            if(frechetSum){
                EE = frechetSE(SE, wr);
            }else{
                EE = expSE(EE);
            }
            At[j] = affine(expSym(SS),Vector3d::Zero()) * EE;
        }else if(blendMode == 2){    //logmatrix3
            Matrix3d G=Matrix3d::Zero();
            Vector3d l=Vector3d::Zero();
            for(int i=0; i<numPrb; i++){
                G +=  wr[i] * logGL[i];
                l += wl[i] * L[i];
            }
            At[j] = ::affine(G.exp(), l);
        }else if(blendMode == 3){   // logmatrix4
            Matrix4d A=Matrix4d::Zero();
            for(int i=0; i<numPrb; i++)
                A +=  wr[i] * logAff[i];
            At[j] = A.exp();
        }else if(blendMode == 5){ // quaternion
            Vector4d q=Vector4d::Zero();
            Matrix3d SS=Matrix3d::Zero();
            Vector3d l=Vector3d::Zero();
            for(int i=0; i<numPrb; i++){
                q += wr[i] * quat[i];
                SS += ws[i] * logS[i];
                l += wl[i] * L[i];
            }
            SS = ::expSym(SS);
            Quaternion<double> Q(q);
            Matrix3d RR = Q.matrix().transpose();
            At[j] = ::affine(SS*RR, l);
        }else if(blendMode==10){
            At[j] = Matrix4d::Zero();
            for(int i=0; i<numPrb; i++){
                At[j] += wr[i] * aff[i];
            }
        }
    }


// compute target vertices position
	MatrixXd G=MatrixXd::Zero(numTet+numPts,3);
    arapG(At, PI, triangles, aff, G);
    MatrixXd Sol = solver.solve(G);
    for(int i=0;i<numPts;i++){
        pts[i].x=Sol(i,0);
        pts[i].y=Sol(i,1);
        pts[i].z=Sol(i,2);
    }
    if(worldMode){
        for(int i=0;i<numPts;i++)
            pts[i] *= localToWorldMatrix.inverse();
    }
    itGeo.setAllPositions(pts);
    return MS::kSuccess;
}

void probeDeformerARAPNode::arapHI(const std::vector<Matrix4d>& PI, const MIntArray& triangles)
{
    int dim = numTet + numPts;
    std::vector<T> tripletList;
    tripletList.reserve(numTet*16);
    Matrix4d Hlist;
	Matrix4d diag=Matrix4d::Identity();
	diag(3,3)=transWeight;
    int s,t;
	for(int i=0;i<numTet;i++){
		Hlist=PI[i].transpose()*diag*PI[i];
		for(int j=0;j<4;j++){
            if(j==3){
                s=numPts+i;
            }else{
                s=triangles[3*i+j];
            }
			for(int k=0;k<4;k++){
                if(k==3){
                    t=numPts+i;
                }else{
                    t=triangles[3*i+k];
                }
                tripletList.push_back(T(t,s,Hlist(j,k)));
			}
		}
	}
    //    // set hard constraint
    //    for(int i=0;i<numConstraint;i++){
    //        tripletList.push_back(T(dim+i,tr[3*constraintTet[i]],1.0/3.0));
    //        tripletList.push_back(T(dim+i,tr[3*constraintTet[i]+1],1.0/3.0));
    //        tripletList.push_back(T(dim+i,tr[3*constraintTet[i]+2],1.0/3.0));
    //        tripletList.push_back(T(tr[3*constraintTet[i]],dim+i,1.0/3.0));
    //        tripletList.push_back(T(tr[3*constraintTet[i]+1],dim+i,1.0/3.0));
    //        tripletList.push_back(T(tr[3*constraintTet[i]+2],dim+i,1.0/3.0));
    //    }
    //    SpMat mat(dim+numConstraint, dim+numConstraint);
    SpMat mat(dim, dim);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    // set soft constraint
    std::vector<T> constraintList;
    constraintList.reserve(numPrb*3);
    F.resize(dim,numPrb);
    F.setZero();
    for(int i=0;i<numPrb;i++){
        constraintList.push_back(T(triangles[3*constraintTet[i]],i,1.0/3.0));
        constraintList.push_back(T(triangles[3*constraintTet[i]+1],i,1.0/3.0));
        constraintList.push_back(T(triangles[3*constraintTet[i]+2],i,1.0/3.0));
    }
    F.setFromTriplets(constraintList.begin(), constraintList.end());
    mat += constraintWeight * F * F.transpose();
    solver.compute(mat);
    assert(solver.info()==Success);
}

void probeDeformerARAPNode::arapG(const std::vector<Matrix4d>& At, const std::vector<Matrix4d>& PI,
                                     const MIntArray& triangles, const std::vector<Matrix4d>& aff, MatrixXd& G)
{
    Matrix4d Glist;
    Matrix4d diag=Matrix4d::Identity();
    diag(3,3)=transWeight;
    for(int i=0;i<numTet;i++){
        Glist=At[i].transpose()*diag*PI[i];
        for(int k=0;k<3;k++){
            for(int j=0;j<3;j++){
                G(triangles[3*i+j],k) += Glist(k,j);
            }
            G(numPts+i,k) += Glist(k,3);
        }
    }
    // set hard constraint
    //	MatrixXf G=MatrixXf::Zero(dim+numConstraint,3);
    //    for(int i=0;i<numConstraint;i++){
    //        RowVector4d cv = constraintVector[i]*aff[i];
    //        G.block(dim+i,0,1,3) << cv(0), cv(1), cv(2);
    //    }
    // set soft constraint
    std::vector<T> constraintList;
    constraintList.reserve(numPrb*3);
    SpMat S(numPrb,3);
    for(int i=0;i<numPrb;i++){
        RowVector4d cv = constraintVector[i]*aff[i];
        constraintList.push_back(T(i,0,cv(0)));
        constraintList.push_back(T(i,1,cv(1)));
        constraintList.push_back(T(i,2,cv(2)));
    }
    S.setFromTriplets(constraintList.begin(), constraintList.end());
    SpMat FS = constraintWeight * F * S;
    G += MatrixXd(FS);
}


// read array of matrix attributes and convert them to Eigen matrices
void probeDeformerARAPNode::readMatrixArray(MArrayDataHandle& handle, std::vector<Matrix4d>& m)
{
    int numPrb=handle.elementCount();
    MMatrix mat;
    for(int i=0;i<numPrb;i++)
    {
        handle.jumpToArrayElement(i);
        mat=handle.inputValue().asMatrix();
        m[i] << mat(0,0), mat(0,1), mat(0,2), mat(0,3),
            mat(1,0), mat(1,1), mat(1,2), mat(1,3),
            mat(2,0), mat(2,1), mat(2,2), mat(2,3),
            mat(3,0), mat(3,1), mat(3,2), mat(3,3);
    }
}

// create nodes
MStatus probeDeformerARAPNode::initialize()
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

    aWeightMode = eAttr.create( "weightMode", "wtm", 0 );
    eAttr.addField( "normal", 0 );
    eAttr.addField( "curve", 1 );
    eAttr.setStorable(true);
    addAttribute( aWeightMode );
    attributeAffects( aWeightMode, outputGeom );

	aWorldMode = nAttr.create( "worldMode", "wrldmd", MFnNumericData::kBoolean, 0 );
    nAttr.setStorable(true);
    addAttribute( aWorldMode );
    attributeAffects( aWorldMode, outputGeom );
    
	aMaxDist = nAttr.create("maxDistance", "md", MFnNumericData::kDouble, 10.0);
    nAttr.setStorable(true);
	addAttribute( aMaxDist );
	attributeAffects( aMaxDist, outputGeom );

	aTransWeight = nAttr.create("translationWeight", "tw", MFnNumericData::kDouble, 0.0001);
    nAttr.setStorable(true);
	addAttribute( aTransWeight );
	attributeAffects( aTransWeight, outputGeom );
    
	aConstraintWeight = nAttr.create("constraintWeight", "cw", MFnNumericData::kDouble, 100.0);
    nAttr.setStorable(true);
	addAttribute( aConstraintWeight );
	attributeAffects( aConstraintWeight, outputGeom );
    
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
MStatus probeDeformerARAPNode::accessoryNodeSetup(MDagModifier& cmd)
{
	MStatus stat;
	MObject thisNode = thisMObject();
	MRampAttribute rWeightCurveR( thisNode, probeDeformerARAPNode::aWeightCurveR, &stat );
	MRampAttribute rWeightCurveS( thisNode, probeDeformerARAPNode::aWeightCurveS, &stat );
	MRampAttribute rWeightCurveL( thisNode, probeDeformerARAPNode::aWeightCurveL, &stat );
	
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




void probeDeformerARAPNode::tetMatrixC(const MPointArray& p, const MIntArray& triangles, std::vector<Matrix4d>& m, std::vector<Vector3d>& tetCenter)
{
    MVector u, v, q;
    for(int i=0;i<numTet;i++)
    {
        tetCenter[i] << (p[triangles[3*i]].x+p[triangles[3*i+1]].x+p[triangles[3*i+2]].x)/3.0,
            (p[triangles[3*i]].y+p[triangles[3*i+1]].y+p[triangles[3*i+2]].y)/3.0,
            (p[triangles[3*i]].z+p[triangles[3*i+1]].z+p[triangles[3*i+2]].z)/3.0;
        u=p[triangles[3*i+1]]-p[triangles[3*i]];
        v=p[triangles[3*i+2]]-p[triangles[3*i]];
        q=u^v;
//        q.normalize();
        
        m[i] << p[triangles[3*i]].x, p[triangles[3*i]].y, p[triangles[3*i]].z, 1,
        p[triangles[3*i+1]].x, p[triangles[3*i+1]].y, p[triangles[3*i+1]].z, 1,
        p[triangles[3*i+2]].x, p[triangles[3*i+2]].y, p[triangles[3*i+2]].z, 1,
        q[0]+p[triangles[3*i]].x,q[1]+p[triangles[3*i]].y, q[2]+p[triangles[3*i]].z,1;
    }
}

// initializer
MStatus initializePlugin( MObject obj )
{
    MStatus status;
    MFnPlugin plugin( obj, "Shizuo KAJI", "0.1", "Any");
    
    status = plugin.registerNode( probeDeformerARAPNode::nodeName, probeDeformerARAPNode::id, probeDeformerARAPNode::creator, probeDeformerARAPNode::initialize, MPxNode::kDeformerNode );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    
    return status;
}

MStatus uninitializePlugin( MObject obj )
{
    MStatus   status;
    MFnPlugin plugin( obj );
    
    status = plugin.deregisterNode( probeDeformerARAPNode::id );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    
    return status;
}

