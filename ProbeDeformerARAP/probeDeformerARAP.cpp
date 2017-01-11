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
MObject probeDeformerARAPNode::aAreaWeighted;
MObject probeDeformerARAPNode::aNeighbourWeighting;

void* probeDeformerARAPNode::creator() { return new probeDeformerARAPNode; }
 
MStatus probeDeformerARAPNode::deform( MDataBlock& data, MItGeometry& itGeo, const MMatrix &localToWorldMatrix, unsigned int mIndex )
{
	MObject thisNode = thisMObject();
    MStatus status;
    MThreadUtils::syncNumOpenMPThreads();    // for OpenMP
    
    bool worldMode = data.inputValue( aWorldMode ).asBool();
    bool areaWeighted = data.inputValue( aAreaWeighted ).asBool();
    short stiffnessMode = data.inputValue( aStiffness ).asShort();
    short blendMode = data.inputValue( aBlendMode ).asShort();
    short tetMode = data.inputValue( aTetMode ).asShort();
    short numIter = data.inputValue( aIteration ).asShort();
    short constraintMode = data.inputValue( aConstraintMode ).asShort();
    short visualisationMode = data.inputValue( aVisualisationMode ).asShort();
    mesh.transWeight = data.inputValue( aTransWeight ).asDouble();
    double constraintWeight = data.inputValue( aConstraintWeight ).asDouble();
    double normExponent = data.inputValue( aNormExponent ).asDouble();
    double visualisationMultiplier = data.inputValue(aVisualisationMultiplier).asDouble();
    MArrayDataHandle hMatrixArray = data.inputArrayValue(aMatrix);
    MArrayDataHandle hInitMatrixArray = data.inputArrayValue(aInitMatrix);
    // check connection
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
    B.setNum(numPrb);
    // read matrices from probes
    std::vector<Matrix4d> initMatrix(numPrb), matrix(numPrb);
    readMatrixArray(hInitMatrixArray, initMatrix);
    readMatrixArray(hMatrixArray, matrix);
    // read vertex positions
    MPointArray Mpts;
    itGeo.allPositions(Mpts);
    int numPts = Mpts.length();
    
    // compute distance
    if(!data.isClean(aARAP) || !data.isClean(aComputeWeight) || isNumProbeChanged){
        // load points list
        if(worldMode){
            for(int j=0; j<numPts; j++ )
                Mpts[j] *= localToWorldMatrix;
        }
        pts.resize(numPts);
        for(int i=0;i<numPts;i++){
            pts[i] << Mpts[i].x, Mpts[i].y, Mpts[i].z;
        }
        // make tetrahedral structure
        getMeshData(data, input, inputGeom, mIndex, tetMode, pts, mesh.tetList, faceList, edgeList, vertexList, mesh.tetMatrix, mesh.tetWeight);
        mesh.dim = removeDegenerate(tetMode, numPts, mesh.tetList, faceList, edgeList, vertexList, mesh.tetMatrix);
        makeTetMatrix(tetMode, pts, mesh.tetList, faceList, edgeList, vertexList, mesh.tetMatrix, mesh.tetWeight);
        makeTetCenterList(tetMode, pts, mesh.tetList, tetCenter);
        mesh.numTet = (int)mesh.tetList.size()/4;
        mesh.computeTetMatrixInverse();
        // initial probe position
        for(int i=0;i<numPrb;i++){
            B.centre[i] = transPart(initMatrix[i]);
        }
        // compute distance between probe and tetrahedra
        D.setNum(numPrb, numPts, mesh.numTet);
        D.computeDistTet(tetCenter, B.centre);
        D.findClosestTet();
        D.computeDistPts(pts, B.centre);
        D.findClosestPts();
        if(!areaWeighted){
            mesh.tetWeight.clear();
            mesh.tetWeight.resize(mesh.numTet,1.0);
        }
    }
    
    // (re)compute ARAP
    if(!data.isClean(aARAP) || isNumProbeChanged){
        // load painted weights
        if(stiffnessMode == SM_PAINT) {
            VectorXd ptsWeight(numPts);
            for (int i=0; !itGeo.isDone(); itGeo.next()){
                double w=weightValue(data, mIndex, itGeo.index());
                ptsWeight[i++] = (w>EPSILON) ? w : EPSILON;
            }
            makeTetWeightList(tetMode, mesh.tetList, faceList, edgeList, vertexList, ptsWeight, mesh.tetWeight);
        }else if(stiffnessMode == SM_LEARN) {
            std::vector<double> tetEnergy(mesh.numTet,0);
            MArrayDataHandle hSupervisedMesh = data.inputArrayValue(aSupervisedMesh);
            int numSupervisedMesh = hSupervisedMesh.elementCount();
            for(int j=0;j<numSupervisedMesh;j++){
                hSupervisedMesh.jumpToElement(j);
                MFnMesh ex_mesh(hSupervisedMesh.inputValue().asMesh());
                MPointArray Mspts;
                ex_mesh.getPoints( Mspts );
                if(numPts != Mspts.length()){
                    MGlobal::displayInfo("incompatible mesh");
                    return MS::kFailure;
                }
                std::vector<Vector3d> spts(numPts);
                for(int i=0;i<numPts;i++){
                    spts[i] << Mspts[i].x, Mspts[i].y, Mspts[i].z;
                }
                std::vector<double> dummy_weight;
                makeTetMatrix(tetMode, spts, mesh.tetList, faceList, edgeList, vertexList, Q, dummy_weight);
                Matrix3d S,R;
                for(int i=0;i<mesh.numTet;i++)  {
                    polarHigham((mesh.tetMatrixInverse[i]*Q[i]).block(0,0,3,3), S, R);
                    tetEnergy[i] += (S-Matrix3d::Identity()).squaredNorm();
                }
            }
            // compute weight (stiffness)
            double max_energy = *std::max_element(tetEnergy.begin(), tetEnergy.end());
            for(int i=0;i<mesh.numTet;i++)  {
                double w = 1.0 - tetEnergy[i]/(max_energy+EPSILON);
                mesh.tetWeight[i] *= w*w;
            }
        }

        // find constraint points
        constraint.resize(3*numPrb);
        for(int i=0;i<numPrb;i++){
            constraint[3*i] = T(i,mesh.tetList[4*D.closestTet[i]],constraintWeight);
            constraint[3*i+1] = T(i,mesh.tetList[4*D.closestTet[i]+1],constraintWeight);
            constraint[3*i+2] = T(i,mesh.tetList[4*D.closestTet[i]+2],constraintWeight);
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
            double constraintRadius = data.inputValue( aConstraintRadius ).asDouble();
            for(int i=0;i<numPrb;i++){
                double r = constraintRadius * probeConstraintRadius[i];
                for(int j=0;j<numPts;j++){
                    if(D.distPts[i][j]<r){
                        constraint.push_back(T(i,j,constraintWeight * pow((r-D.distPts[i][j])/r,normExponent)));
                    }
                }
            }
        }
        int numConstraint=constraint.size();
        mesh.constraintWeight.resize(numConstraint);
        mesh.constraintVal.resize(numConstraint,numPrb);
        for(int cur=0;cur<numConstraint;cur++){
            mesh.constraintWeight[cur] = std::make_pair(constraint[cur].col(), constraint[cur].value());
        }
        //
        isError = mesh.ARAPprecompute();
        status = data.setClean(aARAP);
    }        // END of ARAP precomputation
    
    if(isError>0){
        return MS::kFailure;
    }
    
    // probe weight computation
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
        wr.resize(mesh.numTet);ws.resize(mesh.numTet);wl.resize(mesh.numTet);
        for(int j=0;j<mesh.numTet;j++){
            wr[j].resize(numPrb); ws[j].resize(numPrb); wl[j].resize(numPrb);
        }
        short weightMode = data.inputValue( aWeightMode ).asShort();
        if (weightMode == WM_INV_DISTANCE){
            for(int j=0;j<mesh.numTet;j++){
                double sum=0.0;
                std::vector<double> idist(numPrb);
                for (int i = 0; i<numPrb; i++){
                    idist[i] = probeRadius[i] / pow(D.distTet[i][j], normExponent);
                    sum += idist[i];
                }
                for (int i = 0; i<numPrb; i++){
                    wr[j][i] = ws[j][i] = wl[j][i] = sum > 0 ? idist[i] / sum : 0.0;
                }
            }
        }
        else if (weightMode == WM_CUTOFF_DISTANCE){
            for(int j=0;j<mesh.numTet;j++){
                for (int i = 0; i<numPrb; i++){
                    wr[j][i] = ws[j][i] = wl[j][i] = (D.distTet[i][j] > probeRadius[i])
                    ? 0 : pow((probeRadius[i] - D.distTet[i][j]) / probeRadius[i], normExponent);
                }
            }
        }else if (weightMode == WM_DRAW){
            float val;
            MRampAttribute rWeightCurveR( thisNode, aWeightCurveR, &status );
            MRampAttribute rWeightCurveS( thisNode, aWeightCurveS, &status );
            MRampAttribute rWeightCurveL( thisNode, aWeightCurveL, &status );
            for(int j=0;j<mesh.numTet;j++){
                for (int i = 0; i < numPrb; i++){
                    rWeightCurveR.getValueAtPosition(D.distTet[i][j] / probeRadius[i], val);
                    wr[j][i] = val;
                    rWeightCurveS.getValueAtPosition(D.distTet[i][j] / probeRadius[i], val);
                    ws[j][i] = val;
                    rWeightCurveL.getValueAtPosition(D.distTet[i][j] / probeRadius[i], val);
                    wl[j][i] = val;
                }
            }
        }else if(weightMode & WM_HARMONIC){
            Laplacian harmonicWeighting;
            makeFaceTet(data, input, inputGeom, mIndex, pts, harmonicWeighting.tetList, harmonicWeighting.tetMatrix, harmonicWeighting.tetWeight);
            harmonicWeighting.numTet = (int)harmonicWeighting.tetList.size()/4;
            std::vector<T> weightConstraint(numPrb);
            // the vertex closest to the probe is given probeWeight
            for(int i=0;i<numPrb;i++){
                weightConstraint[i]=T(i,D.closestPts[i],probeWeight[i]);
            }
            // vertices within effectRadius are given probeWeight
            if( data.inputValue( aNeighbourWeighting ).asBool() ){
                for(int i=0;i<numPrb;i++){
                    for(int j=0;j<numPts;j++){
                        if(D.distPts[i][j]<probeRadius[i]){
                            weightConstraint.push_back(T(i,j,probeWeight[i]));
                        }
                    }
                }
            }
            // set boundary condition for weight computation
            int numConstraint=weightConstraint.size();
            harmonicWeighting.constraintWeight.resize(numConstraint);
            harmonicWeighting.constraintVal.resize(numConstraint,numPrb);
            harmonicWeighting.constraintVal.setZero();
            for(int i=0;i<numConstraint;i++){
                harmonicWeighting.constraintVal(i,weightConstraint[i].row())=weightConstraint[i].value();
                harmonicWeighting.constraintWeight[i] = std::make_pair(weightConstraint[i].col(), weightConstraint[i].value());
            }
            // clear tetWeight
            if(!areaWeighted){
                harmonicWeighting.tetWeight.clear();
                harmonicWeighting.tetWeight.resize(harmonicWeighting.numTet,1.0);
            }
            // solve the laplace equation
            if( weightMode == WM_HARMONIC_ARAP){
                harmonicWeighting.computeTetMatrixInverse();
                harmonicWeighting.dim = numPts + harmonicWeighting.numTet;
                isError = harmonicWeighting.ARAPprecompute();
            }else if(weightMode == WM_HARMONIC_COTAN){
                harmonicWeighting.dim = numPts;
                isError = harmonicWeighting.cotanPrecompute();
            }
            if(isError>0) return MS::kFailure;
            std::vector< std::vector<double> > w_tet(numPrb);
            harmonicWeighting.harmonicSolve();
            for(int i=0;i<numPrb;i++){
                makeTetWeightList(tetMode, mesh.tetList, faceList, edgeList, vertexList, harmonicWeighting.Sol.col(i), w_tet[i]);
                for(int j=0;j<mesh.numTet; j++){
                    wr[j][i] = ws[j][i] = wl[j][i] = w_tet[i][j];
                }
            }
        }
        // normalise weights
        short normaliseWeightMode = data.inputValue( aNormaliseWeight ).asShort();
        for(int j=0;j<mesh.numTet;j++){
            D.normaliseWeight(normaliseWeightMode, wr[j]);
            D.normaliseWeight(normaliseWeightMode, ws[j]);
            D.normaliseWeight(normaliseWeightMode, wl[j]);
        }
        status = data.setClean(aComputeWeight);
    } // END of weight computation


    // setting up transformation matrix
    B.rotationConsistency = data.inputValue( aRotationConsistency ).asBool();
    bool frechetSum = data.inputValue( aFrechetSum ).asBool();
    blendedSE.resize(mesh.numTet); blendedR.resize(mesh.numTet); blendedS.resize(mesh.numTet); blendedL.resize(mesh.numTet);A.resize(mesh.numTet);
    for(int i=0;i<numPrb;i++){
        B.Aff[i]=initMatrix[i].inverse()*matrix[i];
    }
    B.parametrise(blendMode);
    

// prepare transform matrix for each simplex
#pragma omp parallel for
	for (int j = 0; j < mesh.numTet; j++){
		// blend matrix
		if (blendMode == BM_SRL){
			blendedS[j] = expSym(blendMat(B.logS, ws[j]));
			Vector3d l = blendMat(B.L, wl[j]);
            blendedR[j] = frechetSum ? frechetSO(B.R, wr[j]) : expSO(blendMat(B.logR, wr[j]));
			A[j] = pad(blendedS[j]*blendedR[j], l);
		}
		else if (blendMode == BM_SSE){
			blendedS[j] = expSym(blendMat(B.logS, ws[j]));
            blendedSE[j] = expSE(blendMat(B.logSE, wr[j]));
			A[j] = pad(blendedS[j], Vector3d::Zero()) * blendedSE[j];
		}
		else if (blendMode == BM_LOG3){
			blendedR[j] = blendMat(B.logGL, wr[j]).exp();
			Vector3d l = blendMat(B.L, wl[j]);
			A[j] = pad(blendedR[j], l);
		}
		else if (blendMode == BM_LOG4){
			A[j] = blendMat(B.logAff, wr[j]).exp();
		}
		else if (blendMode == BM_SQL){
			Vector4d q = blendQuat(B.quat, wr[j]);
			Vector3d l = blendMat(B.L, wl[j]);
			blendedS[j] = blendMatLin(B.S, ws[j]);
			Quaternion<double> RQ(q);
			blendedR[j] = RQ.matrix().transpose();
			A[j] = pad(blendedS[j]*blendedR[j], l);
		}
		else if (blendMode == BM_AFF){
			A[j] = blendMatLin(B.Aff, wr[j]);
		}
	}

    // compute target vertices position
    tetEnergy.resize(mesh.numTet);
    
    // set constraint
    int numConstraints = constraint.size();
    mesh.constraintVal.resize(numConstraints,3);
    RowVector4d cv;
    for(int cur=0;cur<numConstraints;cur++){
        cv = pad(pts[constraint[cur].col()]) * B.Aff[constraint[cur].row()];
        mesh.constraintVal(cur,0) = cv[0];
        mesh.constraintVal(cur,1) = cv[1];
        mesh.constraintVal(cur,2) = cv[2];
    }

    // iterate to determine vertices position
    for(int k=0;k<numIter;k++){
        // solve ARAP
        mesh.ARAPSolve(A);
        // set new vertices position
        new_pts.resize(numPts);
        for(int i=0;i<numPts;i++){
            new_pts[i][0]=mesh.Sol(i,0);
            new_pts[i][1]=mesh.Sol(i,1);
            new_pts[i][2]=mesh.Sol(i,2);
        }
        // if iteration continues
        if(k+1<numIter || visualisationMode == VM_ENERGY){
            std::vector<double> dummy_weight;
            makeTetMatrix(tetMode, new_pts, mesh.tetList, faceList, edgeList, vertexList, Q, dummy_weight);
            Matrix3d S,R,newS,newR;
            if(blendMode == BM_AFF || blendMode == BM_LOG4 || blendMode == BM_LOG3){
                for(int i=0;i<mesh.numTet;i++){
                    polarHigham(A[i].block(0,0,3,3), blendedS[i], blendedR[i]);
                }
            }
            #pragma omp parallel for
            for(int i=0;i<mesh.numTet;i++){
                polarHigham((mesh.tetMatrixInverse[i]*Q[i]).block(0,0,3,3), newS, newR);
                tetEnergy[i] = (newS-blendedS[i]).squaredNorm();
                A[i].block(0,0,3,3) = blendedS[i]*newR;
//                polarHigham((A[i].transpose()*PI[i]*Q[i]).block(0,0,3,3), newS, newR);
//                A[i].block(0,0,3,3) *= newR;
            }
        }
    }
    for(int i=0;i<numPts;i++){
        Mpts[i].x=mesh.Sol(i,0);
        Mpts[i].y=mesh.Sol(i,1);
        Mpts[i].z=mesh.Sol(i,2);
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
            makePtsWeightList(tetMode, numPts, mesh.tetList, faceList, edgeList, vertexList, tetEnergy, ptsColour);
            for(int i=0;i<numPts;i++){
                ptsColour[i] *= visualisationMultiplier;
            }
        }else if(visualisationMode == VM_STIFFNESS){
            makePtsWeightList(tetMode, numPts, mesh.tetList, faceList, edgeList, vertexList, mesh.tetWeight, ptsColour);
            double maxval = *std::max_element(ptsColour.begin(), ptsColour.end());
            for(int i=0;i<numPts;i++){
                ptsColour[i] = 1.0 - ptsColour[i]/maxval;
            }
        }else if(visualisationMode == VM_CONSTRAINT){
            for(int i=0;i<constraint.size();i++){
                ptsColour[constraint[i].col()] += constraint[i].value();
            }
        }else if(visualisationMode == VM_EFFECT){
            std:vector<double> wsum(mesh.numTet);
            for(int j=0;j<mesh.numTet;j++){
                //wsum[j] = std::accumulate(wr[j].begin(), wr[j].end(), 0.0);
                wsum[j]= visualisationMultiplier * wr[j][numPrb-1];
            }
            makePtsWeightList(tetMode, numPts, mesh.tetList, faceList, edgeList, vertexList, wsum, ptsColour);
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
    eAttr.addField( "expSE+expSym", BM_SSE );
    eAttr.addField( "logmatrix3", BM_LOG3 );
    eAttr.addField( "logmatrix4", BM_LOG4 );
    eAttr.addField( "quat+linear", BM_SQL );
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
    
    aNormaliseWeight = eAttr.create( "normaliseWeight", "nw", NM_LINEAR );
    eAttr.addField( "NONE", NM_NONE );
    eAttr.addField( "Linear",  NM_LINEAR );
    eAttr.addField( "Softmax", NM_SOFTMAX );
    eAttr.setStorable(true);
    addAttribute( aNormaliseWeight );
    attributeAffects( aNormaliseWeight, outputGeom );
    attributeAffects( aNormaliseWeight, aComputeWeight );

    aWeightMode = eAttr.create( "weightMode", "wtm", WM_HARMONIC_COTAN );
    eAttr.addField( "inverse", WM_INV_DISTANCE );
    eAttr.addField( "cut-off", WM_CUTOFF_DISTANCE );
    eAttr.addField( "draw", WM_DRAW );
//    eAttr.addField( "harmonic-arap", WM_HARMONIC_ARAP);
    eAttr.addField( "harmonic-cotan", WM_HARMONIC_COTAN);
    eAttr.setStorable(true);
    eAttr.setKeyable(false);
    addAttribute( aWeightMode );
    attributeAffects( aWeightMode, outputGeom );
    attributeAffects( aWeightMode, aComputeWeight );
    
    aConstraintMode = eAttr.create( "constraintMode", "ctm", CONSTRAINT_CLOSEST );
    eAttr.addField( "neighbour",  CONSTRAINT_NEIGHBOUR);
    eAttr.addField( "closestFace", CONSTRAINT_CLOSEST );
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
    
	aTransWeight = nAttr.create("translationWeight", "tw", MFnNumericData::kDouble, 1e-20);
    nAttr.setStorable(true);
	addAttribute( aTransWeight );
	attributeAffects( aTransWeight, outputGeom );
	attributeAffects( aTransWeight, aARAP );

    aAreaWeighted = nAttr.create( "areaWeighted", "aw", MFnNumericData::kBoolean, false );
    nAttr.setStorable(true);
    addAttribute( aAreaWeighted );
    attributeAffects( aAreaWeighted, outputGeom );
    attributeAffects( aAreaWeighted, aARAP );

    aNeighbourWeighting = nAttr.create( "neighbourWeighting", "nghbrw", MFnNumericData::kBoolean, false );
    nAttr.setStorable(true);
    addAttribute( aNeighbourWeighting );
    attributeAffects( aNeighbourWeighting, outputGeom );
    attributeAffects( aNeighbourWeighting, aComputeWeight );
    attributeAffects( aNeighbourWeighting, aARAP );

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

