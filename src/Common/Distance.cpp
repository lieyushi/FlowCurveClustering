#include "Distance.h"


const int& PROCRUSTES_SIZE = 7;

float** distanceMatrix = NULL;

/* ------------------ Compute norm 3 for trajectories ------------------------- */
// given a center trajectory and index of pre-stored vector 
const float getBMetric_3(const VectorXf& row,
						 const int& size,
						 const int& i,
						 const std::vector<std::vector<float> >& rotationSequence
						)
{
	std::vector<float> firstNorm3, secondNorm3;
	getSequence(row, size, firstNorm3);
	secondNorm3 = rotationSequence[i];
	return getBMetric(firstNorm3, secondNorm3);
}

// given two center trajectories for distance measuring 
const float getBMetric_3(const VectorXf& firstRow,
						 const int& size,
						 const VectorXf& secondRow
						)
{
	std::vector<float> firstNorm3, secondNorm3;
	getSequence(firstRow, size, firstNorm3);
	getSequence(secondRow, size, secondNorm3);
	return getBMetric(firstNorm3, secondNorm3);
}


const float getBMetric_3(const int& first,
						 const int& second,
						 const std::vector<std::vector<float> >& rotationSequence
					    )
{
	return getBMetric(rotationSequence[first], rotationSequence[second]);
}

/* -------------------- Finish computing norm 3 for trajectories --------------------*/

/* ------------------ Compute norm 6 for trajectories ------------------------- */
// given a center trajectory and index of pre-stored vector 
const float getBMetric_6(const VectorXf& row,
						 const int& size,
						 const int& i,
						 const std::vector<MultiVariate>& normalMultivariate
						)
{
	MultiVariate centerNormal, neighNormal;
	getNormalMultivariate(row, size, centerNormal);
	neighNormal = normalMultivariate[i];
	return getBMetric(centerNormal, neighNormal);
}

// given two center trajectories for distance measuring 
const float getBMetric_6(const VectorXf& firstRow,
						 const int& size,
						 const VectorXf& secondRow
						)
{
	MultiVariate centerNormal, neighNormal;
	getNormalMultivariate(firstRow, size, centerNormal);
	getNormalMultivariate(secondRow, size, neighNormal);
	return getBMetric(centerNormal, neighNormal);
}

const float getBMetric_6(const int& first,
						 const int& second,
						 const std::vector<MultiVariate>& normalMultivariate
						)
{
	return getBMetric(normalMultivariate[first], normalMultivariate[second]);
}
/* -------------------- Finish computing norm 6 for trajectories --------------------*/


/* ------------------ Compute norm 7 for trajectories ------------------------- */
// given a center trajectory and index of pre-stored vector 
const float getBMetric_7(const VectorXf& row,
						 const int& size,
						 const int& i,
						 const std::vector<std::vector<float> >& rotationSequence
						)
{
	std::vector<float> firstNorm3, secondNorm3;
	getEachFixedSequence(row, size, firstNorm3);
	secondNorm3 = rotationSequence[i];
	return getBMetric(firstNorm3, secondNorm3);
}

// given two center trajectories for distance measuring 
const float getBMetric_7(const VectorXf& firstRow,
						 const int& size,
						 const VectorXf& secondRow
						)
{
	std::vector<float> firstNorm3, secondNorm3;
	getEachFixedSequence(firstRow, size, firstNorm3);
	getEachFixedSequence(secondRow, size, secondNorm3);
	return getBMetric(firstNorm3, secondNorm3);
}

const float getBMetric_7(const int& first,
						 const int& second,
						 const std::vector<std::vector<float> >& rotationSequence
						)
{
	return getBMetric(rotationSequence[first], rotationSequence[second]);
}
/* -------------------- Finish computing norm 7 for trajectories --------------------*/


/* ------------------ Compute norm 9 for trajectories ------------------------- */
// given a center trajectory and index of pre-stored vector 
const float getBMetric_9(const VectorXf& row,
						 const int& size,
						 const int& i,
						 const std::vector<MultiVariate>& normalMultivariate
						)
{
	MultiVariate centerNormal, neighNormal;
	getUnnormalizedMultivariate(row, size, centerNormal);
	neighNormal = normalMultivariate[i];
	return getBMetric(centerNormal, neighNormal);
}

// given two center trajectories for distance measuring 
const float getBMetric_9(const VectorXf& firstRow,
						 const int& size,
						 const VectorXf& secondRow
						)
{
	MultiVariate centerNormal, neighNormal;
	getUnnormalizedMultivariate(firstRow, size, centerNormal);
	getUnnormalizedMultivariate(secondRow, size, neighNormal);
	return getBMetric(centerNormal, neighNormal);
}

const float getBMetric_9(const int& first,
						 const int& second,
						 const std::vector<MultiVariate>& normalMultivariate
						)
{
	return getBMetric(normalMultivariate[first], normalMultivariate[second]);
}

/* -------------------- Finish computing norm 9 for trajectories --------------------*/



const float getBMetric(const std::vector<float>& firstNorm3, 
					   const std::vector<float>& secondNorm3
					  )
{
	float u_a, u_b, sig_a, sig_b, sig_a_inverse, sig_b_inverse, 
		  summation, sum_inverse, tempDist;
	u_a = firstNorm3[0], u_b = secondNorm3[0];
	sig_a = firstNorm3[1], sig_b = secondNorm3[1];
	sig_a *= sig_a, sig_b *= sig_b;
	if(sig_a<=1.0e-8)
	{
		sig_a = 1.0e-8;
		sig_a_inverse = 1.0e8;
	}
	else
		sig_a_inverse = 1.0/sig_a;
	if(sig_b<=1.0e-8)
	{
		sig_b = 1.0e-8;
		sig_b_inverse = 1.0/sig_b;
	}
	summation = sig_a+sig_b;
	sum_inverse = 1.0/summation;
	tempDist = 0.25*log(0.25*(sig_a*sig_b_inverse
			   +sig_b*sig_a_inverse+2))
			   + 0.25*(u_a-u_b)*(u_a-u_b)*sum_inverse;
	return tempDist;
}


const float getBMetric(const MultiVariate& centerNormal, 
					   const MultiVariate& neighNormal
					  ) 
{
	Matrix3f firstCov, secondCov, meanCov, meanCovInverse;
	float sqrtInverse, meanCovDet;
	firstCov = centerNormal.covariance;
	secondCov = neighNormal.covariance;
	meanCov = 0.5*(firstCov+secondCov);
	if(meanCov.determinant()>1.0e-8)
	{
		meanCovInverse = static_cast<Matrix3f>(meanCov.inverse());
		meanCovDet = meanCov.determinant();
	}
	else
	{
		meanCovInverse = pseudoInverse(meanCov);
		meanCovDet = 1.0e8;
	}
	float detMulti = sqrt(firstCov.determinant()*secondCov.determinant());
	sqrtInverse = detMulti>1.0e-8?float(1.0)/detMulti:1.0e8;
	Vector3f meanDiff = centerNormal.meanVec-neighNormal.meanVec;
	float tempDist = 0.125*meanDiff.transpose()*meanCovInverse*meanDiff
			   +0.2*log(meanCovDet*sqrtInverse);
	return tempDist;
}


const float getMetric_10(const VectorXf& centroid,
						 const int& size,
						 const int& index,
						 const std::vector<VectorXf>& unitLength)
{
	const VectorXf& x = unitLength[index];
	VectorXf y(size*3);
	getUnitDirection_byEach(centroid,size,y);
	float length = x.dot(y)/x.size();
	length = min(1.0,(double)length);
	length = max(-1.0,(double)length);
	length = acos(length);
	return length;
}
const float getMetric_10(const VectorXf& firstRow,
						 const int& size,
						 const VectorXf& secondRow)
{
	VectorXf x(size*3), y(size*3);
	getUnitDirection_byEach(firstRow,size,x);
	getUnitDirection_byEach(secondRow,size,y);

	float length = x.dot(y)/x.size();
	length = min(1.0,(double)length);
	length = max(-1.0,(double)length);
	length = acos(length);
	return length;
}
const float getMetric_10(const int& first,
						 const int& second,
						 const std::vector<VectorXf>& unitLength)
{
	const VectorXf& x = unitLength[first];
	const VectorXf& y = unitLength[second];
	float length = x.dot(y)/x.size();
	length = min(1.0,(double)length);
	length = max(-1.0,(double)length);
	length = acos(length);
	return length;
}



const float getNorm(const Eigen::VectorXf& centroid,
					const Eigen::VectorXf& r2,
					const int& index,
					const int& normOption,
					const std::vector<std::vector<float> >& pairwise,
					const std::vector<std::vector<float> >& objectNorm)
{
	assert(centroid.size()==r2.size());
	float length = 0.0;
	switch(normOption)
	{
	case 0:
	default:
		length = (centroid-r2).norm();
		break;

	case 11:
		{
			float dotPro = centroid.dot(r2);
			float firstNorm = centroid.norm();
			float secondNorm = r2.norm();
			float firstInverse, secondInverse;
			if(firstNorm<1.0e-8)
				firstInverse = 1.0e8;
			else
				firstInverse = 1.0/firstNorm;
			if(secondNorm<1.0e-8)
				secondInverse = 1.0e8;
			else
				secondInverse = 1.0/secondNorm;
			dotPro = dotPro*firstInverse*secondInverse;
			dotPro = std::max(dotPro, float(-1.0));
			dotPro = std::min(dotPro, float(1.0));
			length = acos(dotPro)/M_PI;
		}
		break;

	case 1:  /* fraction norm by high-dimensional feature-space */
		{
			for (int i = 0; i < centroid.size(); ++i)
			{
				length += pow(abs(centroid(i)-r2(i)),0.5);
			}
			length = pow(length,2.0);
		}
		break;

	case 2: /* mean value of dot product value, which means it's rotational invariant */
		{
			const int& pointNum = centroid.size()/3-1;
			float dotValue, leftNorm, rightNorm, result;

			std::vector<float> centroidWise;
			std::vector<float> centroidWiseNorm;
			getPairWise_byEach(centroid, pointNum, centroidWise, centroidWiseNorm);

			const std::vector<float>& i_Pairwise = pairwise[index];
			const std::vector<float>& i_PairNorm = objectNorm[index];

			Vector3f left, right;
			for (int i = 0; i < pointNum; ++i)
			{
				leftNorm = centroidWiseNorm[i];
				rightNorm = i_PairNorm[i];
				if(leftNorm >= 1.0e-8 && rightNorm >= 1.0e-8)
				{
					left << centroidWise[3*i], centroidWise[3*i+1], centroidWise[3*i+2];
					right << i_Pairwise[3*i],i_Pairwise[3*i+1],i_Pairwise[3*i+2];
					result = left.dot(right)/*/leftNorm/rightNorm*/;
					result = min(1.0,(double)result);
					result = max(-1.0,(double)result);
					length+=acos(result);
				}
				else
					length+=M_PI/2.0;
			}
			length /= pointNum;
		}
		break;

	case 5: /* rotational invariant line-wise acos angle with normal direction for
			   measuring whether counterclockwise or clockwise orientation */
		{
			const int& pointNum = centroid.size()/3-1;
			float dotValue, leftNorm, rightNorm, normalDot, result;
			Vector3f left, right, normal;

			std::vector<float> centroidWise;
			std::vector<float> centroidNorm;
			getPairWise_byEach(centroid, pointNum, centroidWise, centroidNorm);
			const std::vector<float>& i_Pairwise = pairwise[index];
			const std::vector<float>& i_PairNorm = objectNorm[index];

			left << /*centroid(3)-centroid(0),centroid(4)-centroid(1),centroid(5)-centroid(2)*/
					centroidWise[0], centroidWise[1], centroidWise[2];
			right << /*r2(3)-r2(0),r2(4)-r2(1),r2(5)-r2(2)*/
					i_Pairwise[0], i_Pairwise[1], i_Pairwise[2];
			const Vector3f& Normal = left.cross(right);

			for (int i = 0; i < pointNum; ++i)
			{

				leftNorm = centroidNorm[i];
				rightNorm = i_PairNorm[i];
				if(leftNorm >= 1.0e-8 && rightNorm >= 1.0e-8)
				{
					left << centroidWise[3*i], centroidWise[3*i+1], centroidWise[3*i+2];
					right << i_Pairwise[3*i],i_Pairwise[3*i+1],i_Pairwise[3*i+2];
					result = left.dot(right)/*/leftNorm/rightNorm*/;
					result = min(1.0,(double)result);
					result = max(-1.0,(double)result);
					normal = left.cross(right);
					normalDot = Normal.dot(normal);
					if(normalDot<0)
						length+=-acos(result);
					else
						length+=acos(result);
				}
				else
					length+=M_PI/2.0;
			}
			length /= pointNum;
			length = abs(length);
		}
		break;	

	case 8: /* distance metric defined as mean * standard deviation */
		{
			const int& pointNum = centroid.size()/3-1;
			float dotValue, leftNorm, rightNorm, stdevia = 0.0, angle, result;
			Vector3f left, right;

			std::vector<float> centroidWise;
			std::vector<float> centroidNorm;
			getPairWise_byEach(centroid, pointNum, centroidWise, centroidNorm);

			const std::vector<float>& i_Pairwise = pairwise[index];
			const std::vector<float>& i_PairNorm = objectNorm[index];

			for (int i = 0; i < pointNum; ++i)
			{
				leftNorm = centroidNorm[i];
				rightNorm = i_PairNorm[i];

				if(leftNorm >= 1.0e-8 && rightNorm >= 1.0e-8)
				{
					left << centroidWise[3*i], centroidWise[3*i+1], centroidWise[3*i+2];
					right << i_Pairwise[3*i],i_Pairwise[3*i+1],i_Pairwise[3*i+2];
					result = left.dot(right)/*/leftNorm/rightNorm*/;
					result = min(1.0,(double)result);
					result = max(-1.0,(double)result);
					angle = acos(result);
					length+=angle;
					stdevia+=angle*angle;
				}
				else
				{
					angle=M_PI/2.0;
					length+=angle;
					stdevia+=angle*angle;
				}
			}
			length /= pointNum;
			stdevia = stdevia/pointNum-length*length;
			if(stdevia>0)
				stdevia = sqrt(stdevia/pointNum-length*length);
			else
				stdevia = 1.0e-4;
		}
		break;
	}

	return length;
}


const float getNorm(const VectorXf& centroid,
					const VectorXf& r2,
					const int& firstIndex,
					const int& secondIndex,
					const int& normOption,
					const std::vector<std::vector<float> >& pairwise,
					const std::vector<std::vector<float> >& objectNorm)
{
	assert(centroid.size()==r2.size());
	float length = 0.0;
	switch(normOption)
	{
	case 0:
	default:
		length = (centroid-r2).norm();
		break;

	case 1:  /* fraction norm by high-dimensional feature-space */
		{
			for (int i = 0; i < centroid.size(); ++i)
			{
				length += pow(abs(centroid(i)-r2(i)),0.5);
			}
			length = pow(length,2.0);
		}
		break;

	case 11:
		{
			float dotPro = centroid.dot(r2);
			float firstNorm = centroid.norm();
			float secondNorm = r2.norm();
			float firstInverse, secondInverse;
			if(firstNorm<1.0e-8)
				firstInverse = 1.0e8;
			else
				firstInverse = 1.0/firstNorm;
			if(secondNorm<1.0e-8)
				secondInverse = 1.0e8;
			else
				secondInverse = 1.0/secondNorm;
			dotPro = dotPro*firstInverse*secondInverse;
			dotPro = std::max(dotPro, float(-1.0));
			dotPro = std::min(dotPro, float(1.0));
			length = acos(dotPro)/M_PI;
		}
		break;

	case 2: /* mean value of dot product value, which means it's rotational invariant */
		{
			const int& pointNum = centroid.size()/3-1;
			float dotValue, leftNorm, rightNorm, result;

			const std::vector<float>& i_Pairwise = pairwise[firstIndex];
			const std::vector<float>& j_Pairwise = pairwise[secondIndex];

			const std::vector<float>& i_PairNorm = objectNorm[firstIndex];
			const std::vector<float>& j_PairNorm = objectNorm[secondIndex];

			Vector3f left, right;
			for (int i = 0; i < pointNum; ++i)
			{
				leftNorm = i_PairNorm[i];
				rightNorm = j_PairNorm[i];
				if(leftNorm >= 1.0e-8 && rightNorm >= 1.0e-8)
				{
					left << i_Pairwise[3*i], i_Pairwise[3*i+1], i_Pairwise[3*i+2];
					right << j_Pairwise[3*i],j_Pairwise[3*i+1],j_Pairwise[3*i+2];
					result = left.dot(right)/*/leftNorm/rightNorm*/;
					result = min(1.0,(double)result);
					result = max(-1.0,(double)result);
					length+=acos(result);
				}
				else
					length+=M_PI/2.0;
			}
			length /= pointNum;
		}
		break;

	case 5: /* rotational invariant line-wise acos angle with normal direction for
			   measuring whether counterclockwise or clockwise orientation */
		{
			const int& pointNum = centroid.size()/3-1;
			float dotValue, leftNorm, rightNorm, normalDot, result;
			Vector3f left, right, normal;

			const std::vector<float>& i_Pairwise = pairwise[firstIndex];
			const std::vector<float>& j_Pairwise = pairwise[secondIndex];

			const std::vector<float>& i_PairNorm = objectNorm[firstIndex];
			const std::vector<float>& j_PairNorm = objectNorm[secondIndex];

			left << /*centroid(3)-centroid(0),centroid(4)-centroid(1),centroid(5)-centroid(2)*/
					i_Pairwise[0], i_Pairwise[1], i_Pairwise[2];
			right << /*r2(3)-r2(0),r2(4)-r2(1),r2(5)-r2(2)*/
					j_Pairwise[0], j_Pairwise[1], j_Pairwise[2];
			const Vector3f& Normal = left.cross(right);

			for (int i = 0; i < pointNum; ++i)
			{
				leftNorm = i_PairNorm[i];
				rightNorm = j_PairNorm[i];

				if(leftNorm >= 1.0e-8 && rightNorm >= 1.0e-8)
				{
					left << i_Pairwise[3*i], i_Pairwise[3*i+1], i_Pairwise[3*i+2];
					right << j_Pairwise[3*i],j_Pairwise[3*i+1],j_Pairwise[3*i+2];
					result = left.dot(right)/*/leftNorm/rightNorm*/;
					result = min(1.0,(double)result);
					result = max(-1.0,(double)result);
					normal = left.cross(right);
					normalDot = Normal.dot(normal);
					if(normalDot<0)
						length+=-acos(result);
					else
						length+=acos(result);
				}
				else
					length+=M_PI/2.0;
			}
			length /= pointNum;
			length = abs(length);
		}
		break;	

	case 8: /* distance metric defined as mean * standard deviation */
		{
			const int& pointNum = centroid.size()/3-1;
			float dotValue, leftNorm, rightNorm, stdevia = 0.0, angle, result;
			Vector3f left, right;

			const std::vector<float>& i_Pairwise = pairwise[firstIndex];
			const std::vector<float>& j_Pairwise = pairwise[secondIndex];

			const std::vector<float>& i_PairNorm = objectNorm[firstIndex];
			const std::vector<float>& j_PairNorm = objectNorm[secondIndex];

			for (int i = 0; i < pointNum; ++i)
			{
				leftNorm = i_PairNorm[i];
				rightNorm = j_PairNorm[i];

				if(leftNorm >= 1.0e-8 && rightNorm >= 1.0e-8)
				{
					left << i_Pairwise[3*i], i_Pairwise[3*i+1], i_Pairwise[3*i+2];
					right << j_Pairwise[3*i],j_Pairwise[3*i+1],j_Pairwise[3*i+2];
					result = left.dot(right)/*/leftNorm/rightNorm*/;
					result = min(1.0,(double)result);
					result = max(-1.0,(double)result);
					angle = acos(result);
					length+=angle;
					stdevia+=angle*angle;
				}
				else
				{
					angle=M_PI/2.0;
					length+=angle;
					stdevia+=angle*angle;
				}
			}
			length /= pointNum;
			stdevia = stdevia/pointNum-length*length;
			if(stdevia>0)
				stdevia = sqrt(stdevia/pointNum-length*length);
			else
				stdevia = 1.0e-4;
		}
		break;

	}

	return length;
}


const float getNorm(const Eigen::VectorXf& r1, 
					const Eigen::VectorXf& r2, 
					const int& normOption)
{
	assert(r1.size()==r2.size());
	float length = 0.0;
	switch(normOption)
	{
	case 0:
	default:
		length = (r1-r2).norm();
		break;

	case 1:  /* fraction norm by high-dimensional feature-space */
		{
			for (int i = 0; i < r1.size(); ++i)
			{
				length += pow(abs(r1(i)-r2(i)),0.5);
			}
			length = pow(length,2.0);
		}
		break;

	case 2: /* mean value of dot product value, which means it's rotational invariant */
		{
			const int& pointNum = r1.size()/3-1;
			float dotValue, leftNorm, rightNorm, result;
			Vector3f left, right;
			for (int i = 0; i < pointNum; ++i)
			{
				left << r1(3*i+3)-r1(3*i),r1(3*i+4)-r1(3*i+1),r1(3*i+5)-r1(3*i+2);
				right << r2(3*i+3)-r2(3*i),r2(3*i+4)-r2(3*i+1),r2(3*i+5)-r2(3*i+2);
				dotValue = left.dot(right);
				leftNorm = left.norm(), rightNorm = right.norm();
				if(leftNorm >= 1.0e-8 && rightNorm >=1.0e-8)
				{
					result = dotValue/*/leftNorm/rightNorm*/;
					result = min(1.0,(double)result);
					result = max(-1.0,(double)result);
					length+=acos(result);
				}
				else
					length+=M_PI/2.0;
			}
			length /= pointNum;
		}
		break;

	case 5: /* rotational invariant line-wise acos angle with normal direction for
			   measuring whether counterclockwise or clockwise orientation */
		{
			const int& pointNum = r1.size()/3-1;
			float dotValue, leftNorm, rightNorm, normalDot, result;
			Vector3f left, right, normal;

			left << r1(3)-r1(0),r1(4)-r1(1),r1(5)-r1(2);
			right << r2(3)-r2(0),r2(4)-r2(1),r2(5)-r2(2);

			const Vector3f& Normal = left.cross(right);

			for (int i = 0; i < pointNum; ++i)
			{
				left << r1(3*i+3)-r1(3*i),r1(3*i+4)-r1(3*i+1),r1(3*i+5)-r1(3*i+2);
				right << r2(3*i+3)-r2(3*i),r2(3*i+4)-r2(3*i+1),r2(3*i+5)-r2(3*i+2);
				normal = left.cross(right);
				dotValue = left.dot(right);
				normalDot = Normal.dot(normal);

				leftNorm = left.norm(), rightNorm = right.norm();
				if(leftNorm >= 1.0e-8 && rightNorm >=1.0e-8)
				{
					result = dotValue/*/leftNorm/rightNorm*/;
					result = min(1.0,(double)result);
					result = max(-1.0,(double)result);
					if(normalDot<0)
						length+=-acos(result);
					else
						length+=acos(result);
				}
				else
					length+=M_PI/2.0;
			}
			length /= pointNum;
			length = abs(length);
		}
		break;	

	case 8: /* distance metric defined as mean * standard deviation */
		{
			const int& pointNum = r1.size()/3-1;
			float dotValue, leftNorm, rightNorm, stdevia = 0.0, angle, result;
			Vector3f left, right;
			for (int i = 0; i < pointNum; ++i)
			{
				left << r1(3*i+3)-r1(3*i),r1(3*i+4)-r1(3*i+1),r1(3*i+5)-r1(3*i+2);
				right << r2(3*i+3)-r2(3*i),r2(3*i+4)-r2(3*i+1),r2(3*i+5)-r2(3*i+2);
				dotValue = left.dot(right);
				leftNorm = left.norm(), rightNorm = right.norm();
				if(leftNorm >= 1.0e-8 && rightNorm >=1.0e-8)
				{
					result = dotValue/*/leftNorm/rightNorm*/;
					result = min(1.0,(double)result);
					result = max(-1.0,(double)result);
					angle = acos(result);
					std::cout << angle << std::endl;
					length+=angle;
					stdevia+=angle*angle;
				}
				else
				{
					angle=M_PI/2.0;
					length+=angle;
					stdevia+=angle*angle;
				}
			}
			length /= pointNum;
			stdevia = stdevia/pointNum-length*length;
			if(stdevia>0)
				stdevia = sqrt(stdevia/pointNum-length*length);
			else
				stdevia = 1.0e-4;
			length*=stdevia;
		}
		break;

	case 11:
		{
			float dotPro = r1.dot(r2);
			float firstNorm = r1.norm();
			float secondNorm = r2.norm();
			float firstInverse, secondInverse;
			if(firstNorm<1.0e-8)
				firstInverse = 1.0e8;
			else
				firstInverse = 1.0/firstNorm;
			if(secondNorm<1.0e-8)
				secondInverse = 1.0e8;
			else
				secondInverse = 1.0/secondNorm;
			dotPro = dotPro*firstInverse*secondInverse;
			dotPro = std::max(dotPro, float(-1.0));
			dotPro = std::min(dotPro, float(1.0));
			length = acos(dotPro)/M_PI;
		}
		break;
	}

	return length;
}



const float getDisimilarity(const MatrixXf& data,
							const int& first,
							const int& second,
							const int& normOption,
							const MetricPreparation& object)
{
	return getDisimilarity(data.row(first), data.row(second),
						   first, second, normOption, object);
}


const float getDisimilarity(const VectorXf& others,
							const MatrixXf& data,
							const int& index,
							const int& normOption,
							const MetricPreparation& object)
{
	float length;
	switch(normOption)
	{
	case 0:
	case 1:
	case 2:
	case 5:
	case 8:
	case 11:
		length = getNorm(others, data.row(index),index,normOption,
					object.pairwise, object.pairwiseNorm);
		break;

	case 3:
		length = getBMetric_3(others, others.size()/3-2, index, 
							  object.rotationSequence);
		break;

	case 4:
		length = abs(object.rotation[index]-
				 getRotation(others, others.size()/3-2));
		break;

	case 6:
		length = getBMetric_6(others, others.size()/3-1, index,
							  object.normalMultivariate);
		break;

	case 7:
		length = getBMetric_7(others, others.size()/3-1, index, 
							  object.rotationSequence);
		break;

	case 9:
		length = getBMetric_9(others, others.size()/3-1, index, 
							  object.normalMultivariate);
		break;

	case 10:
		length = getMetric_10(others, others.size()/3, index, 
							  object.unitLength);
		break;

	case 12:
		length = getMetric_MOP(others, data.row(index));
		break;

	case 13:
		length = getMetric_Hausdorff(others, data.row(index));
		break;

	/* signature-based similarity metric with chi-squared test combined with mean-closest */
	case 14:
		length = getSignatureMetric(others,data.row(index),object.pairwise[index]);
		break;

	/* adapted Procrustes distance */
	case 15:
		length = getProcrustesMetric(others, data.row(index));
		//length = getProcrustesMetricSegment(first,second);
		break;

	case 16:
		length = getEntropyMetric(object.pairwise[index], others);
		break;

	default:
		exit(1);
		break;
	}

	return length;
	
}


const float getDisimilarity(const VectorXf& first,
							const VectorXf& second,
							const int& firstIndex,
							const int& secondIndex,
							const int& normOption,
							const MetricPreparation& object)
{
	float length;
	switch(normOption)
	{
	case 0:
	case 1:
	case 2:
	case 5:
	case 8:
	case 11:
		length = getNorm(first, second,firstIndex,secondIndex, normOption,
					object.pairwise, object.pairwiseNorm);
		break;

	case 3:
		length = getBMetric_3(firstIndex, secondIndex,
					object.rotationSequence);
		break;

	case 4:
		length = abs(object.rotation[firstIndex]-
					object.rotation[secondIndex]);
		break;

	case 6:
		length = getBMetric_6(firstIndex, secondIndex, 
					object.normalMultivariate);
		break;

	case 7:
		length = getBMetric_7(firstIndex, secondIndex, 
					object.rotationSequence);
		break;

	case 9:
		length = getBMetric_9(firstIndex, secondIndex,
				    object.normalMultivariate);
		break;

	case 10:
		length = getMetric_10(firstIndex, secondIndex, 
					object.unitLength);
		break;

	case 12:
		length = getMetric_MOP(first, second);
		break;

	case 13:
		length = getMetric_Hausdorff(first, second);
		break;

	case 14:
		length = getSignatureMetric(first,second,object.pairwise[firstIndex],object.pairwise[secondIndex]);
		break;

	case 15:
		length = getProcrustesMetric(first, second);
		//length = getProcrustesMetricSegment(first,second);
		break;

	case 16:
		length = getEntropyMetric(object.pairwise[firstIndex], object.pairwise[secondIndex]);
		break;


	default:
		exit(1);
		break;
	}

	return length;
}


const float getMetric_MOP(const VectorXf& first, const VectorXf& second)
{
	const int& vNum = first.size()/3;
	float result, f_to_s, s_to_f;
	float summation = 0;
	for(int i=0;i<vNum;++i)
	{
		float minDist = FLT_MAX;
		Vector3f m_i = Vector3f(first(3*i),first(3*i+1),first(3*i+2));
		for(int j=0;j<vNum;++j)
		{
			Vector3f n_j = Vector3f(second(3*j),second(3*j+1),second(3*j+2));
			minDist = std::min((m_i-n_j).norm(),minDist);
		}
		summation+=minDist;
	}
	s_to_f = summation/vNum;

	summation = 0;
	for(int i=0;i<vNum;++i)
	{
		float minDist = FLT_MAX;
		Vector3f m_i = Vector3f(second(3*i),second(3*i+1),second(3*i+2));
		for(int j=0;j<vNum;++j)
		{
			Vector3f n_j = Vector3f(first(3*j),first(3*j+1),first(3*j+2));
			minDist = std::min((m_i-n_j).norm(),minDist);
		}
		summation+=minDist;
	}
	f_to_s = summation/vNum;

	result = (f_to_s+s_to_f)/2.0;
	return result;
}

/* get Hausdorff distance between streamlines */
const float getMetric_Hausdorff(const VectorXf& first, const VectorXf& second)
{
	const int& vNum = first.size()/3;
	float result, f_to_s=-1.0, s_to_f=-1.0;
	for(int i=0;i<vNum;++i)
	{
		float minDist = FLT_MAX;
		Vector3f m_i = Vector3f(first(3*i),first(3*i+1),first(3*i+2));
		for(int j=0;j<vNum;++j)
		{
			Vector3f n_j = Vector3f(second(3*j),second(3*j+1),second(3*j+2));
			minDist = std::min((m_i-n_j).norm(),minDist);
		}
		s_to_f=std::max(s_to_f, minDist);
	}

	for(int i=0;i<vNum;++i)
	{
		float minDist = FLT_MAX;
		Vector3f m_i = Vector3f(second(3*i),second(3*i+1),second(3*i+2));
		for(int j=0;j<vNum;++j)
		{
			Vector3f n_j = Vector3f(first(3*j),first(3*j+1),first(3*j+2));
			minDist = std::min((m_i-n_j).norm(),minDist);
		}
		f_to_s=std::max(f_to_s, minDist);
	}

	result = std::max(f_to_s, s_to_f);
	return result;
}



bool getDistanceMatrix(const MatrixXf& data,
				       const int& normOption,
					   const MetricPreparation& object)
{
	try
	{
		const int& Row = data.rows();
		distanceMatrix = new float*[Row];
	#pragma omp parallel for schedule(dynamic) num_threads(8)
		for (int i = 0; i < Row; ++i)
		{
			distanceMatrix[i] = new float[Row];
			for (int j = 0; j < Row; ++j)
			{
				/* don't wish to waste computation on diagonal element */
				if(i==j)
					distanceMatrix[i][j] = 0.0;
				else
					distanceMatrix[i][j] = getDisimilarity(data.row(i), data.row(j), i, j, normOption, object);
			}
		}

		std::cout << "Finished computing distance matrix!" << std::endl;
		return true;
	}
	catch(std::bad_alloc& exc)
	{
		return false;
	}
}


void deleteDistanceMatrix(const int& Row)
{
	if(distanceMatrix)
	{
	#pragma omp parallel for schedule(dynamic) num_threads(8)	
		for (int i = 0; i < Row; ++i)
		{
			if(distanceMatrix[i])
			{
				delete[] distanceMatrix[i];
				distanceMatrix[i] = NULL;
			}	
		}
		delete[] distanceMatrix;
		distanceMatrix = NULL;
	}
}


const float getRotation(const std::vector<vector<float> >& streamline, std::vector<float>& rotation)
{
	if(streamline.empty())
		return -1;
	float result = 0, eachSum;
	const int& size = streamline.size();
	rotation = std::vector<float>(size);
	std::vector<float> eachLine;
	Eigen::Vector3f first, second;
	int lineSize;
	for(int i=0;i<size;++i)
	{
		eachSum = 0;
		eachLine = streamline[i];
		lineSize = eachLine.size()/3-2;
		for(int j=0;j<lineSize;++j)
		{
			first<<eachLine[3*j+3]-eachLine[3*j],eachLine[3*j+4]-eachLine[3*j+1],eachLine[3*j+5]-eachLine[3*j+2];
			second<<eachLine[3*j+6]-eachLine[3*j+3],eachLine[3*j+7]-eachLine[3*j+4],eachLine[3*j+8]-eachLine[3*j+5];
			float angle = first.dot(second)/first.norm()/second.norm();
			angle = std::max(angle,float(-1.0));
			angle = std::min(angle,float(1.0));
			eachSum+=acos(angle);
		}
		rotation[i]=eachSum;
		result+=eachSum;
	}
	result/=size;
	return result;
}


/* get signature-based dissimilarity metric given two elements and their histogram*/
const float getSignatureMetric(const Eigen::VectorXf& firstArray,
							   const Eigen::VectorXf& secondArray,
							   const std::vector<float>& firstHist,
							   const std::vector<float>& secondHist)
{
	/* would choose alpha = 0.5, and 10% of subset vertices for mean_dist */
	const float& Alpha = 0.5;
	const int& SUBSET = 10;

	/* assert whether the size is the same */
	assert(firstArray.size()==secondArray.size());
	assert(firstHist.size()==secondHist.size());

	const int& histSize = firstHist.size();
	const int& vertexCount = firstArray.size()/3;
	const int& size = vertexCount/SUBSET+1;

	Eigen::VectorXf firstSubset(3*size), secondSubset(3*size);


	/* get mean_dist between two sampled subsets */
	int tempPos = 0;
	for(int i=0;i<vertexCount;i+=SUBSET)
	{
		for(int j=0;j<3;++j)
		{
			firstSubset(3*tempPos+j)=firstArray(3*i+j);
			secondSubset(3*tempPos+j)=secondArray(3*i+j);
		}
		++tempPos;
	}

	/* get mean_dist */
	float result = getMetric_MOP(firstSubset, secondSubset);

	float chi_test = 0.0, histDiff, histSum;

	/* get chi_test for two histograms */
	for(int i=0;i<histSize;++i)
	{
		histDiff = firstHist[i]-secondHist[i];
		histSum = firstHist[i]+secondHist[i];
		/* check numerical error */
		if(histSum<1.0e-8)
			continue;

		chi_test+= histDiff*histDiff/histSum;
	}

	/* get combined distance */
	result = (1-Alpha)*chi_test + Alpha*result;

	return result;
}


/* get signature-based dissimilarity metric given centroid */
const float getSignatureMetric(const Eigen::VectorXf& centroid,
							   const Eigen::VectorXf& first,
							   const std::vector<float>& firstHist)
{
	std::vector<float> centroidHist;
	/* get the bin-based histogram for signature */
	getSignatureHist(centroid, BIN_SIZE, centroidHist);

	return getSignatureMetric(centroid,first,centroidHist,firstHist);
}


/* get adapted Procrustes distance */
const float getProcrustesMetric(const Eigen::VectorXf& first,
								const Eigen::VectorXf& second)
{
	assert(first.size()==second.size());

	const int& vertexCount = first.size()/3;

	const int& vertexChanged = vertexCount-2*(PROCRUSTES_SIZE/2);
	const int& newSize = 3*vertexChanged;

	/* assign the segment list */
	Eigen::MatrixXf firstSegment(PROCRUSTES_SIZE,3), secondSegment(PROCRUSTES_SIZE,3), X0;

	int location, rightIndex;

	Eigen::Vector3f first_average, second_average, tempPoint;

	/* A is SVD target, rotation is optimal rotation matrix, and secondPrime is P' after superimposition */
	Eigen::MatrixXf A, rotation, secondPrime = Eigen::MatrixXf(PROCRUSTES_SIZE,3);

	float optimalScaling, traceA, pointDist;

	float result = 0.0;

	/* for all points, assign to them a point set with size of PROCRUSTES_SIZE neighboring points */
	for(int i=0;i<vertexChanged;++i)
	{
		rightIndex = i+PROCRUSTES_SIZE;

		first_average = second_average = Eigen::VectorXf::Zero(3);

		/* get the point set of neighboring 7 points and average */
		for(int j=i;j<rightIndex;++j)
		{
			location = j-i;
			for(int k=0;k<3;++k)
			{
				firstSegment(location,k)=first(3*j+k);
				secondSegment(location,k)=second(3*j+k);
			}

			first_average+=firstSegment.row(location);
			second_average+=secondSegment.row(location);
		}

		first_average/=PROCRUSTES_SIZE;
		second_average/=PROCRUSTES_SIZE;

		/* centralization for the point set */
		for(int j=0;j<PROCRUSTES_SIZE;++j)
		{
			firstSegment.row(j) = firstSegment.row(j)-first_average.transpose();
			secondSegment.row(j) = secondSegment.row(j)-second_average.transpose();
		}

		/* get ssqX and ssqY */
		float ssqX = (firstSegment.cwiseProduct(firstSegment)).sum();
		float ssqY = (secondSegment.cwiseProduct(secondSegment)).sum();

		/* check whether negative or not */
		assert(ssqX > 0 && ssqY > 0);

		ssqX = sqrt(ssqX);
		ssqY = sqrt(ssqY);

		/* reserve the matrix */
		X0 = firstSegment;

		/* scaling for the point set */
		firstSegment/=ssqX;
		secondSegment/=ssqY;

		/* get the optimal rotational matrix by othogonal Procrutes analysis */
		A = firstSegment.transpose()*secondSegment;

		/* perform SVD on A */
		JacobiSVD<MatrixXf> svd(A, ComputeThinU | ComputeThinV);

		/* get the optimal 3D rotation */
		rotation = svd.matrixV()*(svd.matrixU().transpose());

		/* get trace for singular value matrix */
		traceA = svd.singularValues().sum();

		/* get optimal scaling */
		optimalScaling = traceA*ssqX/ssqY;

		/* preset the average to the P' */
		for(int j=0;j<PROCRUSTES_SIZE;++j)
			secondPrime.row(j) = second_average;

		/* get P' in superimposed space */
		secondPrime = ssqX*traceA*secondSegment*A+secondPrime;

		/* compute the distance and store them in the std::vector<float> */
		pointDist = 0.0;
		for(int j=0;j<PROCRUSTES_SIZE;++j)
		{
			tempPoint = X0.row(j)-secondPrime.row(j);
			pointDist+= tempPoint.transpose()*tempPoint;
		}

		/* get the average of P(x,y')^2 */
		result+=pointDist;
	}

	return result/vertexChanged;

}


/* get adapted Procrustes distance. For example, if vec has 100 points, it will calculate mean of 14 points */
const float getProcrustesMetricSegment(const Eigen::VectorXf& first,
									   const Eigen::VectorXf& second)
{
	assert(first.size()==second.size());

	const int& vertexCount = first.size()/3;

	const int& vertexChanged = vertexCount/PROCRUSTES_SIZE;
	const int& newSize = 3*vertexChanged;

	/* assign the segment list */
	Eigen::MatrixXf firstSegment(PROCRUSTES_SIZE,3), secondSegment(PROCRUSTES_SIZE,3), X0;

	int location, rightIndex;

	Eigen::Vector3f first_average, second_average, tempPoint;

	/* A is SVD target, rotation is optimal rotation matrix, and secondPrime is P' after superimposition */
	Eigen::MatrixXf A, rotation, secondPrime = Eigen::MatrixXf(PROCRUSTES_SIZE,3);

	float optimalScaling, traceA, pointDist;

	float result = 0.0;

	/* for all points, assign to them a point set with size of PROCRUSTES_SIZE neighboring points */
	for(int i=0;i<vertexChanged;++i)
	{
		rightIndex = PROCRUSTES_SIZE*i+PROCRUSTES_SIZE;

		first_average = second_average = Eigen::VectorXf::Zero(3);

		/* get the point set of neighboring 7 points and average */
		for(int j=PROCRUSTES_SIZE*i;j<rightIndex;++j)
		{
			location = j-i;
			for(int k=0;k<3;++k)
			{
				firstSegment(location,k)=first(3*j+k);
				secondSegment(location,k)=second(3*j+k);
			}

			first_average+=firstSegment.row(location);
			second_average+=secondSegment.row(location);
		}

		first_average/=PROCRUSTES_SIZE;
		second_average/=PROCRUSTES_SIZE;

		/* centralization for the point set */
		for(int j=0;j<PROCRUSTES_SIZE;++j)
		{
			firstSegment.row(j) = firstSegment.row(j)-first_average.transpose();
			secondSegment.row(j) = secondSegment.row(j)-second_average.transpose();
		}

		/* get ssqX and ssqY */
		float ssqX = (firstSegment.cwiseProduct(firstSegment)).sum();
		float ssqY = (secondSegment.cwiseProduct(secondSegment)).sum();

		/* check whether negative or not */
		assert(ssqX > 0 && ssqY > 0);

		ssqX = sqrt(ssqX);
		ssqY = sqrt(ssqY);

		/* reserve the matrix */
		X0 = firstSegment;

		/* scaling for the point set */
		firstSegment/=ssqX;
		secondSegment/=ssqY;

		/* get the optimal rotational matrix by othogonal Procrutes analysis */
		A = firstSegment.transpose()*secondSegment;

		/* perform SVD on A */
		JacobiSVD<MatrixXf> svd(A, ComputeThinU | ComputeThinV);

		/* get the optimal 3D rotation */
		rotation = svd.matrixV()*(svd.matrixU().transpose());

		/* get trace for singular value matrix */
		traceA = svd.singularValues().sum();

		/* get optimal scaling */
		optimalScaling = traceA*ssqX/ssqY;

		/* preset the average to the P' */
		for(int j=0;j<PROCRUSTES_SIZE;++j)
			secondPrime.row(j) = second_average;

		/* get P' in superimposed space */
		secondPrime = ssqX*traceA*secondSegment*A+secondPrime;

		/* compute the distance and store them in the std::vector<float> */
		pointDist = 0.0;
		for(int j=0;j<PROCRUSTES_SIZE;++j)
		{
			tempPoint = X0.row(j)-secondPrime.row(j);
			pointDist+= tempPoint.transpose()*tempPoint;
		}

		/* get the average of P(x,y')^2 */
		result+=pointDist;
	}

	return result/vertexChanged;
}


/* need to store each label for elements for NID computation */
void generateGroups(const std::vector<std::vector<int> >& storage)
{
	if(storage.empty())
		return;
	std::ofstream readme("../dataset/Storage",ios::out|ios::app);
	if(!readme)
	{
		std::cout << "Error creating Storage!" << std::endl;
		exit(1);
	}

	readme << std::endl;
	const int& groupSize = storage.size();
	std::vector<int> element;
	for(int i=0;i<groupSize;++i)
	{
		element = storage[i];
		if(element.empty())
			continue;
		for(int j=0;j<element.size();++j)
			readme << element[j] << " ";
		readme << std::endl;
	}
	std::cout << std::endl;
	readme.close();
}


/* get illustrative visualization metric for paper An Illustrative Visualization Framework for 3D Vector Fields */
const float getEntropyMetric(const std::vector<float>& firstEntropy,
		                     const std::vector<float>& secondEntropy)
{
	assert(firstEntropy.size()==2);
	assert(secondEntropy.size()==2);

	float first = firstEntropy[0]-secondEntropy[0];
	float second = firstEntropy[1]-secondEntropy[1];

	return sqrt(first*first+second*second);
}

/* get illustrative visualization metric for paper An Illustrative Visualization Framework for 3D Vector Fields */
const float getEntropyMetric(const std::vector<float>& firstEntropy,
		                     const Eigen::VectorXf& array)
{
	assert(firstEntropy.size()==2);

	std::vector<float> secondEntropy;

	getLinearAngularEntropy(array, BUNDLE_SIZE, secondEntropy);

	return getEntropyMetric(firstEntropy, secondEntropy);

}


