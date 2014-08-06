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
 */


#include "StdAfx.h"
#include "probeDeformerARAP.h"

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
#define WM_HARMONIC 3

// constraint mode
#define CONSTRAINT_NEIBOUR 0
#define CONSTRAINT_CLOSEST 1

// stiffness mode
#define SM_NONE 0
#define SM_PAINT 1
#define SM_LEARN 2

// visualisation mode
#define VM_OFF 0
#define VM_ENERGY 1
#define VM_EFFECT 2
#define VM_CONSTRAINT 3
#define VM_STIFFNESS 4

// error codes
#define ERROR_ARAP_PRECOMPUTE 1
#define NUMPRB_AND_ATTR_DIFFERENT 2
#define INCOMPATIBLE_MESH 3

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
MObject probeDeformerARAPNode::aMaxDist;
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
    bool isNumProbeChanged = (numPrb != hMatrixArray.elementCount());
    numPrb = hMatrixArray.elementCount();
    // avoid unnecessary computation
    if(isError>0){
        return MS::kFailure;
    }else if(numPrb != hInitMatrixArray.elementCount() || numPrb == 0 || blendMode == BM_OFF){
        return MS::kSuccess;
    }
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
        status = data.setClean(aARAP);
        // read mesh data
        MArrayDataHandle hInput = data.outputArrayValue( input, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = hInput.jumpToElement( mIndex );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        MObject oInputGeom = hInput.outputValue().child( inputGeom ).asMesh();
        MFnMesh inputGeom(oInputGeom);
        // load points list
        if(worldMode){
            for(int j=0; j<numPts; j++ )
                Mpts[j] *= localToWorldMatrix;
        }
        pts.resize(numPts);
        for(int i=0;i<numPts;i++){
            pts[i] << Mpts[i].x, Mpts[i].y, Mpts[i].z;
        }
        // face list
        MIntArray count, triangles;
        inputGeom.getTriangles( count, triangles );
        int numFaces = triangles.length()/3;
        faceList.resize(3*numFaces);
        for(int i=0;i<3*numFaces;i++){
            faceList[i]=triangles[i];
        }
        // vertex list
        MItMeshPolygon iter(oInputGeom);
        MIntArray faceVertices;
        vertexList.resize(numPts);
        for(int i=0;i<numPts;i++){
            vertexList[i].index = i;
            vertexList[i].connectedTriangles.clear();
        }
        for( ; ! iter.isDone(); iter.next()){
            status = iter.getVertices(faceVertices);
            int count = (int) faceVertices.length();
            for(int j=0;j<count;j++){
                vertexList[faceVertices[j]].connectedTriangles.push_back(faceVertices[(j+1)%count]);
                vertexList[faceVertices[j]].connectedTriangles.push_back(faceVertices[(j+count-1)%count]);
            }
        }
        // set tetrahedra
		std::vector<Matrix4d> P;
        makeEdgeList(faceList, edgeList);
        makeTetList(tetMode, numPts, faceList, edgeList, vertexList, tetList);
        tetMatrix(tetMode, pts, tetList, faceList, edgeList, vertexList, P);
        dim = removeDegenerate(tetMode, numPts, tetList, faceList, edgeList, vertexList, P);
        tetMatrix(tetMode, pts, tetList, faceList, edgeList, vertexList, P);
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
            makeWeightList(tetMode, tetList, faceList, edgeList, vertexList, ptsWeight, tetWeight);
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
                tetMatrix(tetMode, spts, tetList, faceList, edgeList, vertexList, Q);
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
        std::vector< std::vector<double> > d(numPrb);
        closestPts.resize(numPrb);
        constraint.resize(numPrb);
        for(int i=0;i<numPrb;i++){
            constraint[i].clear();
            d[i].resize(numPts);
            closestPts[i] = 0;
            double min_d = HUGE_VAL;
            for(int j=0;j<numPts;j++){
                d[i][j]=(pts[j]-probeCenter[i]).norm();
                if( d[i][j] < min_d){
                    min_d = d[i][j];
                    closestPts[i] = j;
                }
            }
        }
        if( constraintMode == CONSTRAINT_NEIBOUR ){
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
                    if(d[i][j]<r){
                        constraint[i][j] = constraintWeight * pow((r-d[i][j])/r,normExponent);
                    }
                }
            }
        }else if( constraintMode == CONSTRAINT_CLOSEST){
            for(int i=0;i<numPrb;i++){
                constraint[i][closestPts[i]] = constraintWeight;
            }
        }
        

        // prepare ARAP solver
        arapHI(PI, tetList, transWeight);
        if(isError>0) return MS::kFailure;
    }
    
    
    // probe weight computation
    short weightMode = data.inputValue( aWeightMode ).asShort();
    if(!data.isClean(aComputeWeight) || isNumProbeChanged){
        status = data.setClean(aComputeWeight);
        double maxDist = data.inputValue( aMaxDist ).asDouble();
        // load probe weights
        std::vector<double> probeWeight(numPrb), probeRadius(numPrb);
        MArrayDataHandle handle = data.inputArrayValue(aProbeWeight);
        if(handle.elementCount() != numPrb){
            MGlobal::displayInfo("# of Probes and probeWeight are different");
            return MS::kFailure;
        }
        for(int i=0;i<numPrb;i++){
            handle.jumpToArrayElement(i);
            probeWeight[i] = handle.inputValue().asDouble();
            probeRadius[i] = probeWeight[i] * maxDist;
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
                    idist[i] = probeWeight[i] / pow(dist[j][i], normExponent);
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
        }else if(weightMode == WM_HARMONIC){
            harmonicWeight(data,mIndex,probeWeight,tetMode);
            if(isError>0) return MS::kFailure;
        }
        // normalise
        for(int j=0;j<numTet;j++){
            double sum = std::accumulate(wr[j].begin(), wr[j].end(), 0.0);
            if (sum > 1){
                for (int i = 0; i < numPrb; i++){
                    wr[j][i] /= sum;
                }
            }
            sum = std::accumulate(ws[j].begin(), ws[j].end(), 0.0);
            if (sum > 1){
                for (int i = 0; i < numPrb; i++){
                    ws[j][i] /= sum;
                }
            }
            sum = std::accumulate(wl[j].begin(), wl[j].end(), 0.0);
            if (sum > 1){
                for (int i = 0; i < numPrb; i++){
                    wl[j][i] /= sum;
                }
            }
        }
    }


    // read attributes
	bool rotationCosistency = data.inputValue( aRotationConsistency ).asBool();
	bool frechetSum = data.inputValue( aFrechetSum ).asBool();

    // setting transformation matrix
	if( ! rotationCosistency || numPrb != prevNs.size()){
		prevThetas.clear();
		prevThetas.resize(numPrb, 0.0);
		prevNs.clear();
		prevNs.resize(numPrb, Vector3d::Zero());
	}
    std::vector<Matrix3d> logR(numPrb),R(numPrb),logS(numPrb),S(numPrb),logGL(numPrb);
    std::vector<Matrix4d> logSE(numPrb),SE(numPrb),logAff(numPrb),Aff(numPrb);
    std::vector<Vector3d> L(numPrb);
    std::vector<Vector4d> quat(numPrb);
    for(int i=0;i<numPrb;i++)
        Aff[i]=initMatrix[i].inverse()*matrix[i];
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


// prepare transform matrix for each simplex
    std::vector<Matrix4d> At(numTet);

#pragma omp parallel for
	for (int j = 0; j < numTet; j++){
		// blend matrix
		if (blendMode == BM_SRL){
			Matrix3d RR, SS = expSym(blendMat(logS, ws[j]));
			Vector3d l = blendMat(L, wl[j]);
			if (frechetSum){
				RR = frechetSO(R, wr[j]);
			}
			else{
				RR = expSO(blendMat(logR, wr[j]));
			}
			At[j] = pad(SS*RR, l);
		}
		else if (blendMode == BM_SES){
			Matrix4d RR;
			Matrix3d SS = expSym(blendMat(logS, ws[j]));
			if (frechetSum){
				RR = frechetSE(SE, wr[j]);
			}
			else{
				RR = expSE(blendMat(logSE, wr[j]));
			}
			At[j] = pad(SS, Vector3d::Zero()) * RR;
		}
		else if (blendMode == BM_LOG3){
			Matrix3d RR = blendMat(logGL, wr[j]).exp();
			Vector3d l = blendMat(L, wl[j]);
			At[j] = pad(RR, l);
		}
		else if (blendMode == BM_LOG4){
			At[j] = blendMat(logAff, wr[j]).exp();
		}
		else if (blendMode == BM_QSL){
			Vector4d q = blendQuat(quat, wr[j]);
			Vector3d l = blendMat(L, wl[j]);
			Matrix3d SS = blendMatLin(S, ws[j]);
			Quaternion<double> Q(q);
			Matrix3d RR = Q.matrix().transpose();
			At[j] = pad(SS*RR, l);
		}
		else if (blendMode == BM_AFF){
			At[j] = blendMatLin(Aff, wr[j]);
		}
	}


// compute target vertices position
    MatrixXd G(dim,3),Sol;
    std::vector<double> tetEnergy(numTet);
    // iterate to determine vertices position
    for(int k=0;k<numIter;k++){
        // solve ARAP
        arapG(At, PI, tetList, Aff, transWeight, G);
        Sol = solver.solve(G);
        
        // set new vertices position
        std::vector<Vector3d> new_pts(numPts);
        for(int i=0;i<numPts;i++){
            new_pts[i][0]=Sol(i,0);
            new_pts[i][1]=Sol(i,1);
            new_pts[i][2]=Sol(i,2);
        }
        // if iteration continues
        if(k+1<numIter || visualisationMode == VM_ENERGY){
            std::vector<Matrix4d> Q(numTet);
            tetMatrix(tetMode, new_pts, tetList, faceList, edgeList, vertexList, Q);
            Matrix3d S,R,newS,newR;
            for(int i=0;i<numTet;i++){
                polarHigham(At[i].block(0,0,3,3), S, R);
                polarHigham((PI[i]*Q[i]).block(0,0,3,3), newS, newR);
                tetEnergy[i] = (newS-Matrix3d::Identity()).squaredNorm();
                At[i].block(0,0,3,3) = S*newR;
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
                wsum[j] = std::accumulate(wr[j].begin(), wr[j].end(), 0.0);
                wsum[j]= visualisationMultiplier * wr[j][numPrb-1];
            }
            makePtsWeightList(tetMode, numPts, tetList, faceList, edgeList, vertexList, wsum, ptsColour);
        }
        visualise(data, ptsColour);
    }
    
    return MS::kSuccess;
}




/// harmonic weighting
void probeDeformerARAPNode::harmonicWeight(MDataBlock& data, unsigned int mIndex, const std::vector<double>& probeWeight, short tetMode){
    // face list
    MStatus status;
    MArrayDataHandle hInput = data.outputArrayValue( input, &status );
    status = hInput.jumpToElement( mIndex );
    MObject oInputGeom = hInput.outputValue().child( inputGeom ).asMesh();
    MFnMesh inputGeom(oInputGeom);
    MIntArray count, triangles;
    inputGeom.getTriangles( count, triangles );
    std::vector<int> fList(triangles.length());
    for(int i=0;i<triangles.length();i++){
        fList[i]=triangles[i];
    }
    // tet list
    std::vector<int> tList;
    std::vector<Matrix4d> P;
    int numPts=(int)pts.size();
    makeTetList(TM_FACE, numPts, fList, edgeList, vertexList, tList);
    tetMatrix(TM_FACE, pts, tList, fList, edgeList, vertexList, P);
    int num = (int)tList.size()/4;
    // LHS
    std::vector<T> tripletListMat(0);
    tripletListMat.reserve(num*16+2*numPrb);
    Matrix4d Hlist;
	Matrix4d diag=Matrix4d::Identity();
	diag(3,3)=0;
	for(int i=0;i<num;i++){
        P[i]=P[i].inverse().eval();
		Hlist=P[i].transpose()*diag*P[i];
		for(int j=0;j<4;j++){
			for(int k=0;k<4;k++){
                tripletListMat.push_back(T(tList[4*i+j],tList[4*i+k],Hlist(j,k)));
			}
		}
	}
    // set hard constraint
    for(int i=0;i<numPrb;i++){
        tripletListMat.push_back(T(numPts+num+i,closestPts[i],1));
        tripletListMat.push_back(T(closestPts[i],numPts+num+i,1));
    }
    SpMat H(numPts+num+numPrb, numPts+num+numPrb), G(numPts+num+numPrb,numPrb);
    H.setFromTriplets(tripletListMat.begin(), tripletListMat.end());
    // factorise
    SparseLU<SpMat> weightSolver;
//    SpSolver weightSolver;
    //weightSolver.isSymmetric(true);
    weightSolver.compute(H);
    if(weightSolver.info() != Success){
        //        std::string error_mes = solver.lastErrorMessage();
        MGlobal::displayInfo("Cleanup the mesh first: Mesh menu => Cleanup => Remove zero edges, faces");
        isError = ERROR_ARAP_PRECOMPUTE;
    }
    // RHS
    tripletListMat.clear();
    tripletListMat.reserve(numPrb);
    for(int i=0;i<numPrb;i++){
        tripletListMat.push_back(T( numPts+num+i, i, probeWeight[i]));
    }
    G.setFromTriplets(tripletListMat.begin(), tripletListMat.end());
    // solve
    SpMat Sol = weightSolver.solve(G);
    std::vector<double> w_pt(numPts), w_tet(num);
    for (int i=0;i<numPrb; i++){
        for(int j=0;j<numPts; j++){
            w_pt[j] = Sol.coeff(j,i);
        }
        makeWeightList(tetMode, tetList, faceList, edgeList, vertexList, w_pt, w_tet);
        for(int j=0;j<num; j++){
            wr[j][i] = ws[j][i] = wl[j][i] = w_tet[j];
        }
    }
}

///
void probeDeformerARAPNode::arapHI(const std::vector<Matrix4d>& PI, const std::vector<int>& tetList,
                                   double transWeight){
    int numTet = (int)tetList.size()/4;
    std::vector<T> tripletListMat(0);
    tripletListMat.reserve(numTet*16);
    Matrix4d Hlist;
	Matrix4d diag=Matrix4d::Identity();
	diag(3,3)=transWeight;
	for(int i=0;i<numTet;i++){
		Hlist=tetWeight[i] * PI[i].transpose()*diag*PI[i];
		for(int j=0;j<4;j++){
			for(int k=0;k<4;k++){
                tripletListMat.push_back(T(tetList[4*i+j],tetList[4*i+k],Hlist(j,k)));
			}
		}
	}
    SpMat mat(dim, dim);
    mat.setFromTriplets(tripletListMat.begin(), tripletListMat.end());
    // set soft constraint
    int cur=0;
    tripletListMat.resize(0);
    std::vector<T> tripletListFT(0);
    for(int i=0;i<constraint.size();i++){
        std::map<int, double>::iterator iter;
        for(iter = constraint[i].begin(); iter != constraint[i].end(); iter++){
            tripletListMat.push_back(T( iter->first, cur, iter->second));
            tripletListFT.push_back(T( cur, iter->first, 1));
            cur++;
        }
    }
    F.resize(dim,cur);
    F.setFromTriplets(tripletListMat.begin(), tripletListMat.end());
    SpMat FT(cur,dim);
    FT.setFromTriplets(tripletListFT.begin(), tripletListFT.end());
    mat += numTet * F * FT;
    solver.compute(mat);
    if(solver.info() != Success){
//        std::string error_mes = solver.lastErrorMessage();
        MGlobal::displayInfo("Cleanup the mesh first: Mesh menu => Cleanup => Remove zero edges, faces");
        isError = ERROR_ARAP_PRECOMPUTE;
    }
}

void probeDeformerARAPNode::arapG(const std::vector<Matrix4d>& At, const std::vector<Matrix4d>& PI,
                                     const std::vector<int>& tetList,
                                  const std::vector<Matrix4d>& Aff, double transWeight, MatrixXd& G){
    int numTet = (int)tetList.size()/4;
    Matrix4d Glist;
    Matrix4d diag=Matrix4d::Identity();
    diag(3,3)=transWeight;
    G.setZero();
    for(int i=0;i<numTet;i++){
        Glist= tetWeight[i] * PI[i].transpose() * diag * At[i];
        for(int k=0;k<3;k++){
            for(int j=0;j<4;j++){
                G(tetList[4*i+j],k) += Glist(j,k);
            }
        }
    }
    // set soft constraint
    std::vector<T> constraintList(0);
    constraintList.reserve(numPrb*3);
    RowVector4d cv;
    int cur=0;
    for(int i=0;i<numPrb;i++){
        std::map<int, double>::iterator iter;
        for(iter = constraint[i].begin(); iter != constraint[i].end(); iter++){
            cv = pad(pts[iter->first]) * Aff[i];
            constraintList.push_back(T(cur,0,cv[0]));
            constraintList.push_back(T(cur,1,cv[1]));
            constraintList.push_back(T(cur,2,cv[2]));
            cur++;
        }
    }
    SpMat S(cur,3);
    S.setFromTriplets(constraintList.begin(), constraintList.end());
    SpMat FS = numTet * F * S;
    G += MatrixXd(FS);
}

// visualise vertex assigned values
void probeDeformerARAPNode::visualise(MDataBlock& data, std::vector<double>& ptsColour){
    // load target mesh output
    MStatus status;
    MArrayDataHandle outputArray = data.outputArrayValue(outputGeom , &status );
    MDataHandle hOutput = outputArray.inputValue(&status);
    MFnMesh outMesh(hOutput.data());
    MColorArray Colour;
    MIntArray Index;
    // set vertex colour
    for(int i=0;i<ptsColour.size();i++){
        Colour.append( ptsColour[i] ,0,0);
        Index.append(i);
    }
    outMesh.setVertexColors(Colour, Index);
}

// read array of matrix attributes and convert them to Eigen matrices
void probeDeformerARAPNode::readMatrixArray(MArrayDataHandle& handle, std::vector<Matrix4d>& m){
    int numPrb=handle.elementCount();
    m.resize(numPrb);
    MMatrix mat;
    for(int i=0;i<numPrb;i++){
        handle.jumpToArrayElement(i);
        mat=handle.inputValue().asMatrix();
        m[i] << mat(0,0), mat(0,1), mat(0,2), mat(0,3),
            mat(1,0), mat(1,1), mat(1,2), mat(1,3),
            mat(2,0), mat(2,1), mat(2,2), mat(2,3),
            mat(3,0), mat(3,1), mat(3,2), mat(3,3);
    }
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
    addAttribute( aARAP );

    // this attr will be dirtied when weight recomputation is needed
    aComputeWeight = nAttr.create( "computeWeight", "computeWeight", MFnNumericData::kBoolean, true );
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

    aWeightMode = eAttr.create( "weightMode", "wtm", WM_INV_DISTANCE );
    eAttr.addField( "inverse", WM_INV_DISTANCE );
    eAttr.addField( "cut-off", WM_CUTOFF_DISTANCE );
    eAttr.addField( "draw", WM_DRAW );
    eAttr.addField( "harmonic", WM_HARMONIC);
    eAttr.setStorable(true);
    addAttribute( aWeightMode );
    attributeAffects( aWeightMode, outputGeom );
    attributeAffects( aWeightMode, aComputeWeight );
    
    aConstraintMode = eAttr.create( "constraintMode", "ctm", CONSTRAINT_CLOSEST );
    eAttr.addField( "neighbour",  CONSTRAINT_NEIBOUR);
    eAttr.addField( "closestPt", CONSTRAINT_CLOSEST );
    eAttr.setStorable(true);
    addAttribute( aConstraintMode );
    attributeAffects( aConstraintMode, outputGeom );
    attributeAffects( aConstraintMode, aARAP);

    aTetMode = eAttr.create( "tetMode", "tm", TM_FACE );
    eAttr.addField( "face", TM_FACE );
    eAttr.addField( "edge", TM_EDGE );
    eAttr.addField( "vertex", TM_VERTEX );
    eAttr.addField( "vface", TM_VFACE );
    eAttr.setStorable(true);
    addAttribute( aTetMode );
    attributeAffects( aTetMode, outputGeom );
    attributeAffects( aTetMode, aARAP );
    attributeAffects( aTetMode, aComputeWeight );

	aWorldMode = nAttr.create( "worldMode", "wrldmd", MFnNumericData::kBoolean, true );
    nAttr.setStorable(true);
    addAttribute( aWorldMode );
    attributeAffects( aWorldMode, outputGeom );
    attributeAffects( aWorldMode, aARAP );
    
	aMaxDist = nAttr.create("maxDistance", "md", MFnNumericData::kDouble, 8.0);
    nAttr.setStorable(true);
	addAttribute( aMaxDist );
	attributeAffects( aMaxDist, outputGeom );
	attributeAffects( aMaxDist, aComputeWeight );

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

