/**
 * @file ARAP.h
 * @brief ARAP related functions
 * @section LICENSE The MIT License
 * @section requirements:  Eigen library
 * @version 0.10
 * @date  Aug. 2014
 * @author Shizuo KAJI
 */

#pragma once

#include <map>

using namespace Eigen;

typedef SparseMatrix<double> SpMat;
typedef SimplicialLDLT<SpMat> SpSolver;
//    SimplicialLLT<SpMat>;
//    SparseLU<SpMat>;
typedef Triplet<double> T;

#define ERROR_ARAP_PRECOMPUTE 1

///
int ARAPprecompute(const std::vector<Matrix4d>& tetMatrixInverse, const std::vector<int>& tetList,
                    const std::vector<double>& tetWeight, std::vector< std::map<int,double> >& constraint,
                                   double transWeight, int dim, SpMat& constraintMat, SpSolver& solver){
    int numTet = (int)tetList.size()/4;
    std::vector<T> tripletListMat(0);
    tripletListMat.reserve(numTet*16);
    Matrix4d Hlist;
	Matrix4d diag=Matrix4d::Identity();
	diag(3,3)=transWeight;
	for(int i=0;i<numTet;i++){
		Hlist=tetWeight[i] * tetMatrixInverse[i].transpose()*diag*tetMatrixInverse[i];
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
    std::vector<T> tripletListF(0);
    std::map<int, double>::iterator iter;
    for(int i=0;i<constraint.size();i++){
        for(iter = constraint[i].begin(); iter != constraint[i].end(); iter++){
            tripletListMat.push_back(T( iter->first, cur, iter->second));
            tripletListF.push_back(T( cur, iter->first, 1));
            cur++;
        }
    }
    constraintMat.resize(dim,cur);
    constraintMat.setFromTriplets(tripletListMat.begin(), tripletListMat.end());
    SpMat F(cur,dim);
    F.setFromTriplets(tripletListF.begin(), tripletListF.end());
    mat += numTet * constraintMat * F;
    solver.compute(mat);
    if(solver.info() != Success){
        //        std::string error_mes = solver.lastErrorMessage();
        MGlobal::displayInfo("Cleanup the mesh first: Mesh menu => Cleanup => Remove zero edges, faces");
        return ERROR_ARAP_PRECOMPUTE;
    }
    return 0;
}

//
void ARAPSolve(const std::vector<Matrix4d>& A, const std::vector<Matrix4d>& tetMatrixInverse,
                                  const std::vector<int>& tetList,
               const std::vector<double>& tetWeight, const std::vector<Vector3d>& constraintVector,
               double transWeight, int dim, const SpMat& constraintMat, const SpSolver& solver,
               MatrixXd& Sol){
    int numTet = (int)tetList.size()/4;
    Matrix4d Glist;
    Matrix4d diag=Matrix4d::Identity();
    diag(3,3)=transWeight;
    MatrixXd G = MatrixXd::Zero(dim,3);
    for(int i=0;i<numTet;i++){
        Glist= tetWeight[i] * tetMatrixInverse[i].transpose() * diag * A[i];
        for(int k=0;k<3;k++){
            for(int j=0;j<4;j++){
                G(tetList[4*i+j],k) += Glist(j,k);
            }
        }
    }
    // set soft constraint
    std::vector<T> constraintList(0);
    int numConstraint = (int) constraintVector.size();
    constraintList.reserve(numConstraint*3);
    for(int i=0;i<numConstraint;i++){
        constraintList.push_back(T(i,0,constraintVector[i][0]));
        constraintList.push_back(T(i,1,constraintVector[i][1]));
        constraintList.push_back(T(i,2,constraintVector[i][2]));
    }
    SpMat S(numConstraint,3);
    S.setFromTriplets(constraintList.begin(), constraintList.end());
    SpMat FS = numTet * constraintMat * S;
    G += MatrixXd(FS);
    Sol = solver.solve(G);
}
