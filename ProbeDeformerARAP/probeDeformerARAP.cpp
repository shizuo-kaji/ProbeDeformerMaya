/**
 * @file probeDeformerARAP.cpp
 * @brief Probe Deformer ARAP plugin for Maya
 * @section LICENSE The MIT License
 * @section  requirements:  Eigen library, Maya
 * @version 0.15
 * @date  3/Jan/2014
 * @author Shizuo KAJI
 */


#include "StdAfx.h"
#include <set>
#include "probeDeformerARAP.h"
#include "deformerConst.h"

// parametrisation mode
#define BM_SRL 0
#define BM_SES 1
#define BM_LOG3 3
#define BM_LOG4 4
#define BM_QSL 5
#define BM_AFF 10
#define BM_OFF -1


using namespace Eigen;
using namespace AffineLib;
using namespace Tetrise;

MTypeId probeDeformerARAPNode::id( 0x00000104 );
MString probeDeformerARAPNode::nodeName( "probeDeformerARAP" );
MObject probeDeformerARAPNode::aARAP;
MObject probeDeformerARAPNode::aMatrix;
MObject probeDeformerARAPNode::aInitMatrix;
MObject probeDeformerARAPNode::aWorldMode;
MObject probeDeformerARAPNode::aBlendMode;
MObject probeDeformerARAPNode::aTetMode;
MObject probeDeformerARAPNode::aWeightMode;
MObject probeDeformerARAPNode::aWeightCurveR;
MObject probeDeformerARAPNode::aWeightCurveS;
MObject probeDeformerARAPNode::aWeightCurveL;
MObject probeDeformerARAPNode::aEffectRadius;
MObject probeDeformerARAPNode::aTransWeight;
MObject probeDeformerARAPNode::aConstraintWeight;
MObject probeDeformerARAPNode::aRotationConsistency;
MObject probeDeformerARAPNode::aFrechetSum;
MObject probeDeformerARAPNode::aNormExponent;
MObject probeDeformerARAPNode::aIteration;
MObject probeDeformerARAPNode::aConstraintRadius;
MObject probeDeformerARAPNode::aConstraintMode;
MObject probeDeformerARAPNode::aVisualisationMultiplier;
MObject probeDeformerARAPNode::aVisualisationMode;
MObject probeDeformerARAPNode::aSupervisedMesh;
MObject probeDeformerARAPNode::aStiffness;
MObject probeDeformerARAPNode::aProbeWeight;
MObject probeDeformerARAPNode::aProbeConstraintRadius;
MObject probeDeformerARAPNode::aComputeWeight;
MObject probeDeformerARAPNode::aNormaliseWeight;

void* probeDeformerARAPNode::creator() { return new probeDeformerARAPNode; }
 
MStatus probeDeformerARAPNode::deform( MDataBlock& data, MItGeometry& itGeo, const MMatrix &localToWorldMatrix, unsigned int mIndex )
{
	MObject thisNode = thisMObject();
    MStatus status;
    MThreadUtils::syncNumOpenMPThreads();    // for OpenMP
    
    bool worldMode = data.inputValue( aWorldMode ).asBool();
    short stiffnessMode = data.inputValue( aStiffness ).asShort();
    short blendMode = data.inputValue( aBlendMode ).asShort();
    short tetMode = data.inputValue( aTetMode ).asShort();
    short numIter = data.inputValue( aIteration ).asShort();
    short constraintMode = data.inputValue( aConstraintMode ).asShort();
    short visualisationMode = data.inputValue( aVisualisationMode ).asShort();
    double transWeight = data.inputValue( aTransWeight ).asDouble();
    double constraintWeight = data.inputValue( aConstraintWeight ).asDouble();
    double normExponent = data.inputValue( aNormExponent ).asDouble();
    double constraintRadius = data.inputValue( aConstraintRadius ).asDouble();
    double visualisationMultiplier = data.inputValue(aVisualisationMultiplier).asDouble();
    MArrayDataHandle hMatrixArray = data.inputArrayValue(aMatrix);
    MArrayDataHandle hInitMatrixArray = data.inputArrayValue(aInitMatrix);
    // avoid unnecessary computation
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
        deleteAttr(data, aProbeConstraintRadius, indices);
        deleteAttr(data, aProbeWeight, indices);
    }
    bool isNumProbeChanged = (numPrb != hMatrixArray.elementCount());
    numPrb = hMatrixArray.elementCount();
    // read matrices
    std::vector<Matrix4d> initMatrix(numPrb), matrix(numPrb);
    readMatrixArray(hInitMatrixArray, initMatrix);
    readMatrixArray(hMatrixArray, matrix);
    MPointArray Mpts;
    itGeo.allPositions(Mpts);
    int numPts = Mpts.length();
    int numTet = (int)tetList.size()/4;
    // (re)compute ARAP
    if(!data.isClean(aARAP) || isNumProbeChanged){
        // load points list
        if(worldMode){
            for(int j=0; j<numPts; j++ )
                Mpts[j] *= localToWorldMatrix;
        }
        pts.resize(numPts);
        for(int i=0;i<numPts;i++){
            pts[i] << Mpts[i].x, Mpts[i].y, Mpts[i].z;
        }
        // set tetrahedra
		std::vector<Matrix4d> P;
        getMeshData(data, input, inputGeom, mIndex, tetMode, pts, tetList, faceList, edgeList, vertexList, P);
        dim = removeDegenerate(tetMode, numPts, tetList, faceList, edgeList, vertexList, P);
        makeTetMatrix(tetMode, pts, tetList, faceList, edgeList, vertexList, P);
        makeTetCenterList(tetMode, pts, tetList, tetCenter);
        numTet = (int)tetList.size()/4;
        PI.resize(numTet);
		for(int i=0;i<numTet;i++){
			PI[i] = P[i].inverse().eval();
        }

        // load painted weights
        if(stiffnessMode == SM_PAINT) {
            std::vector<double> ptsWeight(numPts);
            for (int i=0; !itGeo.isDone(); itGeo.next()){
                double w=weightValue(data, mIndex, itGeo.index());
                ptsWeight[i++] = (w>EPSILON) ? w : EPSILON;
            }
            makeTetWeightList(tetMode, tetList, faceList, edgeList, vertexList, ptsWeight, tetWeight);
        }else if(stiffnessMode == SM_LEARN) {
            tetWeight.resize(numTet);
            std::vector<double> tetEnergy(numTet,0);
            MArrayDataHandle hSupervisedMesh = data.inputArrayValue(aSupervisedMesh);
            int numSupervisedMesh = hSupervisedMesh.elementCount();
            for(int j=0;j<numSupervisedMesh;j++){
                hSupervisedMesh.jumpToElement(j);
                MFnMesh mesh(hSupervisedMesh.inputValue().asMesh());
                MPointArray Mspts;
                mesh.getPoints( Mspts );
                if(numPts != Mspts.length()){
                    MGlobal::displayInfo("incompatible mesh");
                    return MS::kFailure;
                }
                std::vector<Vector3d> spts(numPts);
                for(int i=0;i<numPts;i++){
                    spts[i] << Mspts[i].x, Mspts[i].y, Mspts[i].z;
                }
                std::vector<Matrix4d> Q(numTet);
                makeTetMatrix(tetMode, spts, tetList, faceList, edgeList, vertexList, Q);
                Matrix3d S,R;
                for(int i=0;i<numTet;i++)  {
                    polarHigham((PI[i]*Q[i]).block(0,0,3,3), S, R);
                    tetEnergy[i] += (S-Matrix3d::Identity()).squaredNorm();
                }
            }
            // compute weight (stiffness)
            double max_energy = *std::max_element(tetEnergy.begin(), tetEnergy.end());
            for(int i=0;i<numTet;i++)  {
                double w = 1.0 - tetEnergy[i]/(max_energy+EPSILON);
                tetWeight[i] = w*w;
            }
        }else{
            tetWeight.clear();
            tetWeight.resize(numTet,1.0);
        }

        // probe position
        probeCenter.resize(numPrb);
        for(int i=0;i<numPrb;i++)
            probeCenter[i] << initMatrix[i](3,0),initMatrix[i](3,1),initMatrix[i](3,2);
        
        // compute distance between probe and tetrahedra
        dist.resize(numTet);
        for(int j=0;j<numTet;j++){
            dist[j].resize(numPrb);
            for(int i=0; i<numPrb; i++){
                dist[j][i] = (tetCenter[j]-probeCenter[i]).norm();
            }
        }

        // find constraint points
        distPts.resize(numPrb);
        closestPts.resize(numPrb);
        constraint.resize(numPrb);
        for(int i=0;i<numPrb;i++){
            constraint[i].clear();
            distPts[i].resize(numPts);
            closestPts[i] = 0;
            double min_d = HUGE_VAL;
            for(int j=0;j<numPts;j++){
                distPts[i][j]=(pts[j]-probeCenter[i]).norm();
                if( distPts[i][j] < min_d){
                    min_d = distPts[i][j];
                    closestPts[i] = j;
                }
            }
        }
        if( constraintMode == CONSTRAINT_NEIGHBOUR ){
            std::vector<double> probeConstraintRadius(numPrb);
            MArrayDataHandle handle = data.inputArrayValue(aProbeConstraintRadius);
            if(handle.elementCount() != numPrb){
                MGlobal::displayInfo("# of Probes and probeConstraintRadius are different");
                return MS::kFailure;
            }
            for(int i=0;i<numPrb;i++){
                handle.jumpToArrayElement(i);
                probeConstraintRadius[i]=handle.inputValue().asDouble();
            }
            for(int i=0;i<numPrb;i++){
                double r = constraintRadius * probeConstraintRadius[i];
                for(int j=0;j<numPts;j++){
                    if(distPts[i][j]<r){
                        constraint[i][j] = constraintWeight * pow((r-distPts[i][j])/r,normExponent);
                    }
                }
            }
        }else if( constraintMode == CONSTRAINT_CLOSEST){
            for(int i=0;i<numPrb;i++){
                constraint[i][closestPts[i]] = constraintWeight;
            }
        }
        // prepare ARAP solver
        isError = ARAPprecompute(PI, tetList, tetWeight, constraint, transWeight, dim, constraintMat, solver);
        // END of precomputation
        status = data.setClean(aARAP);
    }
    if(isError>0){
        return MS::kFailure;
    }
    
    // probe weight computation
    short weightMode = data.inputValue( aWeightMode ).asShort();
    if(!data.isClean(aComputeWeight) || isNumProbeChanged){
        // load probe weights
        MArrayDataHandle handle = data.inputArrayValue(aProbeWeight);
        if(handle.elementCount() != numPrb){
            MGlobal::displayInfo("# of Probes and probeWeight are different");
            isError = ERROR_ATTR;
            return MS::kFailure;
        }
        double effectRadius = data.inputValue( aEffectRadius ).asDouble();
        std::vector<double> probeWeight(numPrb), probeRadius(numPrb);
        for(int i=0;i<numPrb;i++){
            handle.jumpToArrayElement(i);
            probeWeight[i] = handle.inputValue().asDouble();
            probeRadius[i] = probeWeight[i] * effectRadius;
        }
        wr.resize(numTet);ws.resize(numTet);wl.resize(numTet);
        for(int j=0;j<numTet;j++){
            wr[j].resize(numPrb); ws[j].resize(numPrb); wl[j].resize(numPrb);
        }
        if (weightMode == WM_INV_DISTANCE){
            for(int j=0;j<numTet;j++){
                double sum=0.0;
                std::vector<double> idist(numPrb);
                for (int i = 0; i<numPrb; i++){
                    idist[i] = probeRadius[i] / pow(dist[j][i], normExponent);
                    sum += idist[i];
                }
                assert(sum > 0);
                for (int i = 0; i<numPrb; i++){
                    wr[j][i] = ws[j][i] = wl[j][i] = idist[i] / sum;
                }
            }
        }
        else if (weightMode == WM_CUTOFF_DISTANCE){
            for(int j=0;j<numTet;j++){
                for (int i = 0; i<numPrb; i++){
                    wr[j][i] = ws[j][i] = wl[j][i] = (dist[j][i] > probeRadius[i])
                    ? 0 : pow((probeRadius[i] - dist[j][i]) / probeRadius[i], normExponent);
                }
            }
        }else if (weightMode == WM_DRAW){
            float val;
            MRampAttribute rWeightCurveR( thisNode, aWeightCurveR, &status );
            MRampAttribute rWeightCurveS( thisNode, aWeightCurveS, &status );
            MRampAttribute rWeightCurveL( thisNode, aWeightCurveL, &status );
            for(int j=0;j<numTet;j++){
                for (int i = 0; i < numPrb; i++){
                    rWeightCurveR.getValueAtPosition(dist[j][i] / probeRadius[i], val);
                    wr[j][i] = val;
                    rWeightCurveS.getValueAtPosition(dist[j][i] / probeRadius[i], val);
                    ws[j][i] = val;
                    rWeightCurveL.getValueAtPosition(dist[j][i] / probeRadius[i], val);
                    wl[j][i] = val;
                }
            }
        }else if(weightMode == WM_HARMONIC || weightMode == WM_HARMONIC_NEIBOUR){
            std::vector<int> fList,tList;
            std::vector< std::vector<double> > ptsWeight(numPrb), w_tet(numPrb);
            std::vector<Matrix4d> P;
            int d=makeFaceTet(data, input, inputGeom, mIndex, pts, fList, tList, P);
            std::vector< std::map<int,double> > weightConstraint(numPrb);
            std::vector<double> weightConstraintValue(0);
            for(int i=0;i<numPrb;i++){
                weightConstraint[i].clear();
            }
            if( weightMode == WM_HARMONIC_NEIBOUR ){
                for(int i=0;i<numPrb;i++){
                    for(int j=0;j<numPts;j++){
                        if(distPts[i][j]<effectRadius){
                            weightConstraint[i][j] = 1;
                            weightConstraintValue.push_back(probeWeight[i]);
                        }
                    }
                }
            }else if( weightMode == WM_HARMONIC){
                for(int i=0;i<numPrb;i++){
                    weightConstraint[i][closestPts[i]] = 1;
                    weightConstraintValue.push_back(probeWeight[i]);
                }
            }
            isError = harmonicWeight(d, P, tList, fList, weightConstraint, weightConstraintValue, ptsWeight);
            if(isError>0) return MS::kFailure;
            for(int i=0;i<numPrb;i++){
                makeTetWeightList(tetMode, tetList, faceList, edgeList, vertexList, ptsWeight[i], w_tet[i]);
                for(int j=0;j<numTet; j++){
                    wr[j][i] = ws[j][i] = wl[j][i] = w_tet[i][j];
                }
            }
        }
        // normalise
        bool normaliseWeight = data.inputValue( aNormaliseWeight ).asBool();
        for(int j=0;j<numTet;j++){
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


    // read attributes
	bool rotationCosistency = data.inputValue( aRotationConsistency ).asBool();
	bool frechetSum = data.inputValue( aFrechetSum ).asBool();

    // setting transformation matrix
	if( ! rotationCosistency || numPrb != logSE.size() || numPrb != logR.size()){
		logSE.clear();
		logSE.resize(numPrb, Matrix4d::Zero().eval());
		logR.clear();
		logR.resize(numPrb, Matrix3d::Zero().eval());
    }
    
    SE.resize(numPrb);
    logAff.resize(numPrb); Aff.resize(numPrb);
    R.resize(numPrb); logS.resize(numPrb); S.resize(numPrb); logGL.resize(numPrb);
    L.resize(numPrb); quat.resize(numPrb); A.resize(numTet); blendedSE.resize(numTet);
    blendedR.resize(numTet); blendedS.resize(numTet); blendedL.resize(numTet);
    for(int i=0;i<numPrb;i++)
        Aff[i]=initMatrix[i].inverse()*matrix[i];
    if(blendMode == BM_SRL || blendMode == BM_SES || blendMode == BM_QSL){
        for(int i=0;i<numPrb;i++){
            parametriseGL(Aff[i].block(0,0,3,3), logS[i] ,R[i]);
            L[i] = transPart(Aff[i]);
            if(blendMode == BM_SRL){
                logR[i]=logSOc(R[i], logR[i]);
                S[i]=expSym(logS[i]);
            }else if(blendMode == BM_SES){
                SE[i]=pad(R[i], L[i]);
                logSE[i]=logSEc(SE[i], logSE[i]);
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


// prepare transform matrix for each simplex
#pragma omp parallel for
	for (int j = 0; j < numTet; j++){
		// blend matrix
		if (blendMode == BM_SRL){
			blendedS[j] = frechetSum ? frechetSym(S, ws[j]) : expSym(blendMat(logS, ws[j]));
			Vector3d l = blendMat(L, wl[j]);
            blendedR[j] = frechetSum ? frechetSO(R, wr[j]) : expSO(blendMat(logR, wr[j]));
			A[j] = pad(blendedS[j]*blendedR[j], l);
		}
		else if (blendMode == BM_SES){
			blendedS[j] = expSym(blendMat(logS, ws[j]));
            blendedSE[j] = expSE(blendMat(logSE, wr[j]));
			A[j] = pad(blendedS[j], Vector3d::Zero()) * blendedSE[j];
		}
		else if (blendMode == BM_LOG3){
			blendedR[j] = blendMat(logGL, wr[j]).exp();
			Vector3d l = blendMat(L, wl[j]);
			A[j] = pad(blendedR[j], l);
		}
		else if (blendMode == BM_LOG4){
			A[j] = blendMat(logAff, wr[j]).exp();
		}
		else if (blendMode == BM_QSL){
			Vector4d q = blendQuat(quat, wr[j]);
			Vector3d l = blendMat(L, wl[j]);
			blendedS[j] = blendMatLin(S, ws[j]);
			Quaternion<double> Q(q);
			blendedR[j] = Q.matrix().transpose();
			A[j] = pad(blendedS[j]*blendedR[j], l);
		}
		else if (blendMode == BM_AFF){
			A[j] = blendMatLin(Aff, wr[j]);
		}
	}


    // compute target vertices position
    MatrixXd Sol;
    tetEnergy.resize(numTet);
    
    // set constraint
    std::vector<Vector3d> constraintVector(0);
    constraintVector.reserve(numPrb * numPts);
    RowVector4d cv;
    Vector3d cvv;
    for(int i=0;i<numPrb;i++){
        std::map<int, double>::iterator iter;
        for(iter = constraint[i].begin(); iter != constraint[i].end(); iter++){
            cv = pad(pts[iter->first]) * Aff[i];
            cvv << cv[0], cv[1], cv[2];
            constraintVector.push_back(cvv);
        }
    }
    
    // iterate to determine vertices position
    for(int k=0;k<numIter;k++){
        // solve ARAP
        ARAPSolve(A, PI, tetList, tetWeight, constraintVector, transWeight, dim, constraintMat, solver, Sol);
        // set new vertices position
        new_pts.resize(numPts);
        for(int i=0;i<numPts;i++){
            new_pts[i][0]=Sol(i,0);
            new_pts[i][1]=Sol(i,1);
            new_pts[i][2]=Sol(i,2);
        }
        // if iteration continues
        if(k+1<numIter || visualisationMode == VM_ENERGY){
            Q.resize(numTet);
            makeTetMatrix(tetMode, new_pts, tetList, faceList, edgeList, vertexList, Q);
            Matrix3d S,R,newS,newR;
            if(blendMode == BM_AFF || blendMode == BM_LOG4 || blendMode == BM_LOG3){
                for(int i=0;i<numTet;i++){
                    polarHigham(A[i].block(0,0,3,3), blendedS[i], blendedR[i]);
                }
            }
            #pragma omp parallel for
            for(int i=0;i<numTet;i++){
                polarHigham((PI[i]*Q[i]).block(0,0,3,3), newS, newR);
                tetEnergy[i] = (newS-blendedS[i]).squaredNorm();
                A[i].block(0,0,3,3) = blendedS[i]*newR;
//                polarHigham((A[i].transpose()*PI[i]*Q[i]).block(0,0,3,3), newS, newR);
//                A[i].block(0,0,3,3) *= newR;
            }
        }
    }
    for(int i=0;i<numPts;i++){
        Mpts[i].x=Sol(i,0);
        Mpts[i].y=Sol(i,1);
        Mpts[i].z=Sol(i,2);
    }
    if(worldMode){
        for(int i=0;i<numPts;i++)
            Mpts[i] *= localToWorldMatrix.inverse();
    }
    itGeo.setAllPositions(Mpts);
    
    // set vertex colour
    if(visualisationMode != VM_OFF){
        std::vector<double> ptsColour(numPts, 0.0);
        if(visualisationMode == VM_ENERGY){
            makePtsWeightList(tetMode, numPts, tetList, faceList, edgeList, vertexList, tetEnergy, ptsColour);
            for(int i=0;i<numPts;i++){
                ptsColour[i] *= visualisationMultiplier;
            }
        }else if(visualisationMode == VM_STIFFNESS){
            makePtsWeightList(tetMode, numPts, tetList, faceList, edgeList, vertexList, tetWeight, ptsColour);
            double maxval = *std::max_element(ptsColour.begin(), ptsColour.end());
            for(int i=0;i<numPts;i++){
                ptsColour[i] = 1.0 - ptsColour[i]/maxval;
            }
        }else if(visualisationMode == VM_CONSTRAINT){
            for(int i=0;i<constraint.size();i++){
                std::map<int, double>::iterator iter;
                for(iter = constraint[i].begin(); iter != constraint[i].end(); iter++){
                    ptsColour[iter->first] += iter->second;
                }
            }
        }else if(visualisationMode == VM_EFFECT){
            std:vector<double> wsum(numTet);
            for(int j=0;j<numTet;j++){
                //wsum[j] = std::accumulate(wr[j].begin(), wr[j].end(), 0.0);
                wsum[j]= visualisationMultiplier * wr[j][numPrb-1];
            }
            makePtsWeightList(tetMode, numPts, tetList, faceList, edgeList, vertexList, wsum, ptsColour);
        }
        visualise(data, outputGeom, ptsColour);
    }
    
    return MS::kSuccess;
}




// create attributes
MStatus probeDeformerARAPNode::initialize(){
    MFnTypedAttribute tAttr;
    MFnNumericAttribute nAttr;
    MFnEnumAttribute eAttr;
    MFnMatrixAttribute mAttr;
   	MRampAttribute rAttr;

    // this attr will be dirtied when ARAP recomputation is needed
    aARAP = nAttr.create( "arap", "arap", MFnNumericData::kBoolean, true );
    nAttr.setStorable(false);
    nAttr.setKeyable(false);
    nAttr.setHidden(true);
    addAttribute( aARAP );

    // this attr will be dirtied when weight recomputation is needed
    aComputeWeight = nAttr.create( "computeWeight", "computeWeight", MFnNumericData::kBoolean, true );
    nAttr.setStorable(false);
    nAttr.setKeyable(false);
    nAttr.setHidden(true);
    addAttribute( aComputeWeight );

    aMatrix = mAttr.create("probeMatrix", "pm");
    mAttr.setStorable(false);
    mAttr.setHidden(true);
    mAttr.setArray(true);
    mAttr.setUsesArrayDataBuilder(true);
    mAttr.setDisconnectBehavior(MFnMatrixAttribute::kDelete);
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
    eAttr.setKeyable(false);
    addAttribute( aBlendMode );
    attributeAffects( aBlendMode, outputGeom );

	aRotationConsistency = nAttr.create( "rotationConsistency", "rc", MFnNumericData::kBoolean, false );
    nAttr.setKeyable(false);
    nAttr.setStorable(true);
    addAttribute( aRotationConsistency );
    attributeAffects( aRotationConsistency, outputGeom );

	aFrechetSum = nAttr.create( "frechetSum", "fs", MFnNumericData::kBoolean, false );
    nAttr.setKeyable(false);
    nAttr.setStorable(true);
    addAttribute( aFrechetSum );
    attributeAffects( aFrechetSum, outputGeom );
    
    aNormaliseWeight = nAttr.create( "normaliseWeight", "nw", MFnNumericData::kBoolean, true );
    nAttr.setStorable(true);
    addAttribute( aNormaliseWeight );
    attributeAffects( aNormaliseWeight, outputGeom );
    attributeAffects( aNormaliseWeight, aComputeWeight );

    aWeightMode = eAttr.create( "weightMode", "wtm", WM_HARMONIC );
    eAttr.addField( "inverse", WM_INV_DISTANCE );
    eAttr.addField( "cut-off", WM_CUTOFF_DISTANCE );
    eAttr.addField( "draw", WM_DRAW );
    eAttr.addField( "harmonic-closest", WM_HARMONIC);
    eAttr.addField( "harmonic-neibour", WM_HARMONIC_NEIBOUR);
    eAttr.setStorable(true);
    eAttr.setKeyable(false);
    addAttribute( aWeightMode );
    attributeAffects( aWeightMode, outputGeom );
    attributeAffects( aWeightMode, aComputeWeight );
    
    aConstraintMode = eAttr.create( "constraintMode", "ctm", CONSTRAINT_CLOSEST );
    eAttr.addField( "neighbour",  CONSTRAINT_NEIGHBOUR);
    eAttr.addField( "closestPt", CONSTRAINT_CLOSEST );
    eAttr.setStorable(true);
    eAttr.setKeyable(false);
    addAttribute( aConstraintMode );
    attributeAffects( aConstraintMode, outputGeom );
    attributeAffects( aConstraintMode, aARAP);

    aTetMode = eAttr.create( "tetMode", "tm", TM_FACE );
    eAttr.addField( "face", TM_FACE );
    eAttr.addField( "edge", TM_EDGE );
    eAttr.addField( "vertex", TM_VERTEX );
    eAttr.addField( "vface", TM_VFACE );
    eAttr.setStorable(true);
    eAttr.setKeyable(false);
    addAttribute( aTetMode );
    attributeAffects( aTetMode, outputGeom );
    attributeAffects( aTetMode, aARAP );
    attributeAffects( aTetMode, aComputeWeight );

	aWorldMode = nAttr.create( "worldMode", "wrldmd", MFnNumericData::kBoolean, true );
    nAttr.setStorable(true);
    nAttr.setKeyable(false);
    addAttribute( aWorldMode );
    attributeAffects( aWorldMode, outputGeom );
    attributeAffects( aWorldMode, aARAP );
    
	aEffectRadius = nAttr.create("effectRadius", "er", MFnNumericData::kDouble, 8.0);
    nAttr.setMin( EPSILON );
    nAttr.setStorable(true);
	addAttribute( aEffectRadius );
	attributeAffects( aEffectRadius, outputGeom );
	attributeAffects( aEffectRadius, aComputeWeight );

	aTransWeight = nAttr.create("translationWeight", "tw", MFnNumericData::kDouble, 0.0001);
    nAttr.setStorable(true);
	addAttribute( aTransWeight );
	attributeAffects( aTransWeight, outputGeom );
	attributeAffects( aTransWeight, aARAP );
    
	aConstraintWeight = nAttr.create("constraintWeight", "cw", MFnNumericData::kDouble, 1.0);
    nAttr.setStorable(true);
	addAttribute( aConstraintWeight );
	attributeAffects( aConstraintWeight, outputGeom );
	attributeAffects( aConstraintWeight, aARAP );
    
	aNormExponent = nAttr.create("normExponent", "ne", MFnNumericData::kDouble, 1.0);
    nAttr.setStorable(true);
	addAttribute( aNormExponent );
	attributeAffects( aNormExponent, outputGeom );
	attributeAffects( aNormExponent, aARAP );
	attributeAffects( aNormExponent, aComputeWeight );
    
	aIteration = nAttr.create("iteration", "it", MFnNumericData::kShort, 1);
    nAttr.setStorable(true);
    addAttribute(aIteration);
    attributeAffects(aIteration, outputGeom);
    
	aConstraintRadius = nAttr.create("constraintRadius", "cr", MFnNumericData::kDouble, 1.0);
    nAttr.setStorable(true);
	addAttribute( aConstraintRadius );
	attributeAffects( aConstraintRadius, outputGeom );
	attributeAffects( aConstraintRadius, aARAP );

    aVisualisationMode = eAttr.create( "visualisationMode", "vm", VM_OFF );
    eAttr.addField( "off", VM_OFF );
    eAttr.addField( "energy", VM_ENERGY );
    eAttr.addField( "effect", VM_EFFECT );
    eAttr.addField( "constraint", VM_CONSTRAINT );
    eAttr.addField( "stiffness", VM_STIFFNESS );
    eAttr.setStorable(true);
    addAttribute( aVisualisationMode );
    attributeAffects( aVisualisationMode, outputGeom );
    
	aVisualisationMultiplier = nAttr.create("visualisationMultiplier", "vmp", MFnNumericData::kDouble, 1.0);
    nAttr.setStorable(true);
	addAttribute( aVisualisationMultiplier );
	attributeAffects( aVisualisationMultiplier, outputGeom );
    
    aStiffness = eAttr.create( "stiffnessMode", "stfm", SM_NONE );
    eAttr.addField( "off", SM_NONE );
    eAttr.addField( "painted weight", SM_PAINT );
    eAttr.addField( "learn", SM_LEARN );
    eAttr.setStorable(true);
    addAttribute( aStiffness );
    attributeAffects( aStiffness, outputGeom );
    attributeAffects( aStiffness, aARAP );
    
	aSupervisedMesh = tAttr.create("supervisedMesh", "svmesh", MFnData::kMesh);
    tAttr.setStorable(true);
    tAttr.setArray(true);
    tAttr.setUsesArrayDataBuilder(true);
    addAttribute(aSupervisedMesh);
	attributeAffects( aSupervisedMesh, outputGeom );
	attributeAffects( aSupervisedMesh, aARAP );
    
    aProbeWeight = nAttr.create("probeWeight", "prw", MFnNumericData::kDouble, 1.0);
    nAttr.setArray(true);
    nAttr.setStorable(true);
    nAttr.setUsesArrayDataBuilder(true);
    addAttribute(aProbeWeight);
	attributeAffects( aProbeWeight, outputGeom );
	attributeAffects( aProbeWeight, aComputeWeight );

    aProbeConstraintRadius = nAttr.create("probeConstraintRadius", "prcr", MFnNumericData::kDouble, 1.0);
    nAttr.setArray(true);
    nAttr.setStorable(true);
    nAttr.setUsesArrayDataBuilder(true);
    addAttribute(aProbeConstraintRadius);
	attributeAffects( aProbeConstraintRadius, outputGeom );
	attributeAffects( aProbeConstraintRadius, aARAP );

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
    
    // Make the deformer weights paintable
    MGlobal::executeCommand( "makePaintable -attrType multiFloat -sm deformer probeDeformerARAP weights;" );

    return MS::kSuccess;
}

// create ramp attributes
MStatus probeDeformerARAPNode::accessoryNodeSetup(MDagModifier& cmd){
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


// this deformer also changes colours
void probeDeformerARAPNode::postConstructor(){
	setDeformationDetails(kDeformsColors);
}


// initializer
MStatus initializePlugin( MObject obj ){
    MStatus status;
    MFnPlugin plugin( obj, "Shizuo KAJI", "0.1", "CREST");
    
    status = plugin.registerNode( probeDeformerARAPNode::nodeName, probeDeformerARAPNode::id, probeDeformerARAPNode::creator, probeDeformerARAPNode::initialize, MPxNode::kDeformerNode );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    
    return status;
}

MStatus uninitializePlugin( MObject obj ){
    MStatus   status;
    MFnPlugin plugin( obj );
    
    status = plugin.deregisterNode( probeDeformerARAPNode::id );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    
    return status;
}

