/*
 * evrot.cpp
 *  https://github.com/pthimon/clustering
 *  Created on: 04-Mar-2009
 *      Author: sbutler
 */

#include "Evrot.h"

#include <map>

Evrot::Evrot(const Eigen::MatrixXf& X, int method):
	mMethod(method),
	mNumDims(X.cols()),
	mNumData(X.rows()),
	mNumAngles((int)(mNumDims*(mNumDims-1)/2)), // get the number of angles
	ik(Eigen::VectorXi(mNumAngles)),
	jk(Eigen::VectorXi(mNumAngles)),
	mX(X),
	mClusters(std::vector<std::vector<int> >(mNumDims)) //allocate clusters vector
{
	// build index mapping (to index upper triangle)
	int k = 0;
	for( int i=0; i<mNumDims-1; i++ ){
		for( int j=i+1; j<=mNumDims-1; j++ ){
			ik[k] = i;
			jk[k] = j;
			k++;
		}
	}

	evrot();
}



Evrot::~Evrot()
{

}



void Evrot::evrot()
{

	// definitions
	int max_iter = 100;
	float dQ,Q,Q_new,Q_old1,Q_old2,Q_up,Q_down;
	float alpha;
	int iter,d;

	Eigen::VectorXf theta = Eigen::VectorXf::Zero(mNumAngles);
	Eigen::VectorXf theta_new = Eigen::VectorXf::Zero(mNumAngles);

	Q = evqual(mX); // initial quality

	Q_old1 = Q;
	Q_old2 = Q;
	iter = 0;

	while( iter < max_iter ){ // iterate to refine quality
		iter++;
		for( d = 0; d < mNumAngles; d++ ){
			if( mMethod == 2 ){ // descend through numerical drivative
				alpha = 0.1;
				{
					// move up
					theta_new[d] = theta[d] + alpha;
					Eigen::MatrixXf Xrot = rotate_givens(theta_new);
					Q_up = evqual(Xrot);
				}
				{
					// move down
					theta_new[d] = theta[d] - alpha;
					Eigen::MatrixXf Xrot = rotate_givens(theta_new);
					Q_down = evqual(Xrot);
				}

				// update only if at least one of them is better
				if( Q_up > Q || Q_down > Q){
					if( Q_up > Q_down ){
						theta[d] = theta[d] + alpha;
						theta_new[d] = theta[d];
						Q = Q_up;
					} else {
						theta[d] = theta[d] - alpha;
						theta_new[d] = theta[d];
						Q = Q_down;
					}
				}
			} else { // descend through true derivative
				alpha = 1.0;
				dQ = evqualitygrad(theta, d);
				theta_new[d] = theta[d] - alpha * dQ;
				Eigen::MatrixXf Xrot = rotate_givens(theta_new);
				Q_new = evqual(Xrot);
				if( Q_new > Q){
					theta[d] = theta_new[d];
					Q = Q_new;
				}
				else{
					theta_new[d] = theta[d];
				}
			}
		}
		// stopping criteria
		if( iter > 2 ){
			if( Q - Q_old2 < 1e-3 ){
				break;
			}
		}
		Q_old2 = Q_old1;
		Q_old1 = Q;
	}


	mXrot = rotate_givens(theta_new);
	cluster_assign();

	//output
	mQuality = Q;
}



void Evrot::cluster_assign()
{
	// find max of each row
	Eigen::VectorXi max_index_col(mNumData);
#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int i=0; i<mNumData; i++ )
	{
		int col=0;
		float mValue = mXrot.row(i).cwiseAbs().maxCoeff(&col);
		/*
		int row, col;
		mXrot.row(i).cwise().abs().maxCoeff(&row, &col);
		*/
		max_index_col[i] = col;
	}

	// prepare cluster assignments
#pragma omp parallel for schedule(dynamic) num_threads(8)
	for(int j=0; j<mNumDims; j++ )
	{  // loop over all columns
		for(int i=0; i<mNumData; i++ )
		{ // loop over all rows
			if( max_index_col[i] == j ){
				mClusters[j].push_back(i);
			}
		}
	}

/* delete cluster that has zero elements in case that vanishing vector won't create trouble */
	std::vector<std::vector<int> > tempCluster;
	for(int i=0;i<mClusters.size();++i)
		if(!mClusters[i].empty())
			tempCluster.push_back(mClusters[i]);
	mClusters.clear();
	mClusters = tempCluster;
}


float Evrot::evqual(const Eigen::MatrixXf& X)
{
	// take the square of all entries and find max of each row
	Eigen::MatrixXf X2(X.rows(), X.cols());
#pragma omp parallel for schedule(dynamic) num_threads(8)
	for(int i=0;i<X.rows();++i)
	{
		for(int j=0;j<X.cols();++j)
		{
			X2(i,j)=X(i,j)*X(i,j);
		}
	}

	Eigen::VectorXf max_values(X.rows());

#pragma omp parallel for schedule(dynamic) num_threads(8)
	for(int i=0;i<X.rows();++i)
		max_values(i)=X2.row(i).maxCoeff();

	// compute cost
#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int i=0; i<mNumData; i++ )
	{
		X2.row(i) = X2.row(i) / max_values[i];
	}

	float J = 1.0 - (X2.sum()/mNumData -1.0)/mNumDims;

	return J;
}



float Evrot::evqualitygrad(const Eigen::VectorXf& theta, const int& angle_index)
{
	// build V,U,A
	Eigen::MatrixXf V = gradU(theta, angle_index);

	Eigen::MatrixXf U1 = build_Uab(theta, 0,angle_index-1);
	Eigen::MatrixXf U2 = build_Uab(theta, angle_index+1,mNumAngles-1);

	Eigen::MatrixXf A = mX*U1*V*U2;

	// rotate vecs according to current angles
	Eigen::MatrixXf Y = rotate_givens(theta);

	// find max of each row
	Eigen::VectorXf max_values(mNumData);
	Eigen::VectorXi max_index_col(mNumData);
#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int i=0; i<mNumData; i++ ) {
		int row, col;
		Y.row(i).cwiseAbs().maxCoeff(&row, &col);
		max_values[i] = Y(i,col);
		max_index_col[i] = col;
	}

	// compute gradient
	float dJ=0, tmp1, tmp2;
	for( int j=0; j<mNumDims; j++ ){  // loop over all columns
		for( int i=0; i<mNumData; i++ ){ // loop over all rows
			tmp1 = A(i,j) * Y(i,j) / (max_values[i]*max_values[i]);
			tmp2 = A(i,max_index_col[i]) * (Y(i,j)*Y(i,j)) / (max_values[i]*max_values[i]*max_values[i]);
			dJ += tmp1-tmp2;
		}
	}
	dJ = 2*dJ/mNumData/mNumDims;

	return dJ;
}

Eigen::MatrixXf Evrot::rotate_givens(const Eigen::VectorXf& theta)
{
	Eigen::MatrixXf G = build_Uab(theta, 0, mNumAngles-1);
	Eigen::MatrixXf Y = mX*G;
	return Y;
}

Eigen::MatrixXf Evrot::build_Uab(const Eigen::VectorXf& theta, const int& a, const int& b)
{
	int k,i;
	//set Uab to be an identity matrix
	Eigen::MatrixXf Uab(mNumDims,mNumDims);
	Uab.setZero();
	Uab.setIdentity();

	if( b < a ) {
		return Uab;
	}

	float tt,u_ik;
	for( k=a; k<=b; k++ ){
		tt = theta[k];
	#pragma omp parallel for schedule(dynamic) num_threads(8)
		for( i=0; i<mNumDims; i++ ) {
			u_ik = 			Uab(i,ik[k]) * cos(tt) - Uab(i,jk[k]) * sin(tt);
			Uab(i,jk[k]) = 	Uab(i,ik[k]) * sin(tt) + Uab(i,jk[k]) * cos(tt);
			Uab(i,ik[k]) = u_ik;
		}
	}
	return Uab;
}

Eigen::MatrixXf Evrot::gradU(const Eigen::VectorXf& theta, const int& k)
{
	Eigen::MatrixXf V(mNumDims,mNumDims);
	V.setZero();

	V(ik[k],ik[k]) = -sin(theta[k]);
	V(ik[k],jk[k]) = cos(theta[k]);
	V(jk[k],ik[k]) = -cos(theta[k]);
	V(jk[k],jk[k]) = -sin(theta[k]);

	return V;
}
