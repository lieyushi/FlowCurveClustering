#include "PreComputing.h"

void getSequence(const VectorXf& array, 
				 const int& size, 
				 std::vector<float>& rowSequence)
{
	rowSequence = std::vector<float>(2);
	float dotValue, leftNorm, rightNorm, meanRotation = 0.0, deviation = 0.0, angle, result;
	Vector3f left, right;
	for (int j = 0; j < size; ++j)
	{
		left << array(j*3+3)-array(j*3), array(j*3+4)-array(j*3+1), array(j*3+5)-array(j*3+2);
		right << array(j*3+6)-array(j*3+3), array(j*3+7)-array(j*3+4), array(j*3+8)-array(j*3+5);
		dotValue = left.dot(right);
		leftNorm = left.norm();
		rightNorm = right.norm();
		if(leftNorm >= 1.0e-8 && rightNorm >=1.0e-8)
		{
			result = dotValue/leftNorm/rightNorm;
			result = min(1.0,(double)result);
			result = max(-1.0,(double)result);
			angle = acos(result);
			meanRotation += angle;
			deviation += angle*angle;
		}
	}
	meanRotation /= size;
	rowSequence[0] = meanRotation;
	int stdDevia = deviation/size-(meanRotation*meanRotation);
	if(stdDevia<0)
		stdDevia = 1.0e-8;
	rowSequence[1] = sqrt(stdDevia);
}

const float getRotation(const VectorXf& array, 
						const int& size)
{
	float dotValue, leftNorm, rightNorm, meanRotation = 0.0, result;
	Vector3f left, right;
	for (int j = 0; j < size; ++j)
	{
		left << array(j*3+3)-array(j*3), array(j*3+4)-array(j*3+1), array(j*3+5)-array(j*3+2);
		right << array(j*3+6)-array(j*3+3), array(j*3+7)-array(j*3+4), array(j*3+8)-array(j*3+5);
		dotValue = left.dot(right);
		leftNorm = left.norm();
		rightNorm = right.norm();
		if(leftNorm >= 1.0e-8 && rightNorm >=1.0e-8)
		{
			result = dotValue/leftNorm/rightNorm;
			result = min(1.0,(double)result);
			result = max(-1.0,(double)result);
			meanRotation += acos(result);
		}
	}
	meanRotation/=size;
	return meanRotation;
}


void getNormalMultivariate(const VectorXf& array, 
				 	 	   const int& size, 
				 	 	   MultiVariate& rowSequence)
{
	MatrixXf normalDirection(size,3);
	float leftNorm;
	Vector3f left;
	VectorXf unitOne(size);
	for (int j = 0; j < size; ++j)
	{
		left << array(j*3+3)-array(j*3), array(j*3+4)-array(j*3+1), array(j*3+5)-array(j*3+2);
		leftNorm = left.norm();
		if(leftNorm >= 1.0e-8)
		{
			for(int k=0;k<3;k++)
				/* record each line segment normal direction */
				normalDirection(j,k) = left(k)/leftNorm;
		}
		else
		{
			for(int k=0;k<3;k++)
				/* if norm is small, mark them as zero to tell identical points */
				normalDirection(j,k) = 0.0;
		}
		unitOne(j) = 1.0;
	}

	VectorXf meanNormal(3);
	for (int i = 0; i < 3; ++i)
	{
		meanNormal(i) = normalDirection.transpose().row(i).mean();
	}

	MatrixXf tempMatrix = normalDirection-unitOne*meanNormal.transpose();
	rowSequence.covariance = tempMatrix.transpose()*tempMatrix/(size-1);
	rowSequence.meanVec = meanNormal;
}



void getEachFixedSequence(const VectorXf& array, 
				 		  const int& size, 
				 		  std::vector<float>& rowSequence)
{
	rowSequence = std::vector<float>(2);
	float dotValue, leftNorm, meanRotation = 0.0, deviation = 0.0, angle, result;
	Vector3f left, xRay;
	xRay << 1.0,0.0,0.0;
	for (int j = 0; j < size; ++j)
	{
		left << array(j*3+3)-array(j*3), array(j*3+4)-array(j*3+1), array(j*3+5)-array(j*3+2);
		dotValue = left.dot(xRay);
		leftNorm = left.norm();
		if(leftNorm >= 1.0e-8)
		{
			result = dotValue/leftNorm;
			result = min(1.0,(double)result);
			result = max(-1.0,(double)result);
			angle = acos(result);
			meanRotation += angle;
			deviation += angle*angle;
		}
		else
		{
			angle = M_PI/2.0;
			meanRotation += angle;
			deviation += angle*angle;
		}
	}
	meanRotation /= size;
	rowSequence[0] = meanRotation;
	int stdDevia = deviation/size-(meanRotation*meanRotation);
	if(stdDevia<0)
		stdDevia = 1.0e-8;
	rowSequence[1] = sqrt(stdDevia);
}


void getUnnormalizedMultivariate(const VectorXf& array, 
				 	 	  		 const int& size, 
				 	 	  		 MultiVariate& rowSequence)
{
	MatrixXf normalDirection(size,3);
	float leftNorm;
	Vector3f left;
	VectorXf unitOne(size);
	for (int j = 0; j < size; ++j)
	{
		left << array(j*3+3)-array(j*3), array(j*3+4)-array(j*3+1), array(j*3+5)-array(j*3+2);
		leftNorm = left.norm();
		if(leftNorm >= 1.0e-8)
		{
			for(int k=0;k<3;k++)
				/* record each line segment normal direction */
				normalDirection(j,k) = left(k);
		}
		else
		{
			for(int k=0;k<3;k++)
				/* if norm is small, mark them as zero to tell identical points */
				normalDirection(j,k) = 0.0;
		}
		unitOne(j) = 1.0;
	}

	VectorXf meanNormal(3);
	for (int i = 0; i < 3; ++i)
	{
		meanNormal(i) = normalDirection.transpose().row(i).mean();
	}

	MatrixXf tempMatrix = normalDirection-unitOne*meanNormal.transpose();
	rowSequence.covariance = tempMatrix.transpose()*tempMatrix/(size-1);
	rowSequence.meanVec = meanNormal;
}


void getUnitDirection_byEach(const VectorXf& array, 
							 const int& pointNum, 
							 VectorXf& direction)
{
	Vector3f left;
	float leftNorm;
	for (int i = 0; i < pointNum; ++i)
	{
		left << array(3*i), array(3*i+1), array(3*i+2);
		leftNorm = left.norm();
		// I Know it's hardly possible to have smaller norm, but just in case
		if(leftNorm>=1.0e-8)
		{
			for (int j = 0; j < 3; ++j)
			{
				direction(3*i+j) = left(j)/leftNorm;
			}
		}
		else
		{
			for (int j = 0; j < 3; ++j)
			{
				direction(3*i+j) = 0;
			}
		}
	}
}							 


void getPairWise_byEach(const VectorXf& data,
						const int& size,
					 	std::vector<float>& wiseVec,
					 	std::vector<float>& wiseNorm)
{
	if(wiseVec.empty())
		wiseVec = std::vector<float>(3*size);
	
	if(wiseNorm.empty())
		wiseNorm = std::vector<float>(size);
	
	for (int i = 0; i < size; ++i)
	{
		float leftNorm;
		Vector3f left;
		left << data(3*i+3)-data(3*i),data(3*i+4)-data(3*i+1),data(3*i+5)-data(3*i+2);
		leftNorm = left.norm();
		if(leftNorm >= 1.0e-8)
		{
			for (int j = 0; j < 3; ++j)
			{
				wiseVec[3*i+j] = left(j)/leftNorm;
			}
			wiseNorm[i] = leftNorm;
		}
		else
		{
			for (int j = 0; j < 3; ++j)
			{
				wiseVec[3*i+j] = 0.0;
			}
			wiseNorm[i] = 0.0;
		}
	}
}


/* get the bin-based histogram for signature */
void getSignatureHist(const Eigen::VectorXf& array,
					  const int& binNum,
					  std::vector<float>& histogram)
{
	/* if empty vector, should allocate memory ahead of time */
	if(histogram.empty())
		histogram = std::vector<float>(binNum);

	/* get how many vertices you'll have */
	const int& segmentNum = array.size()/3-1;

	/* how many vertices on each bin on average */
	const int& binSize = segmentNum/binNum;

	/* first several has binSize+1 vertices, while the rest have binSize vertices */
	const int& residueNum = segmentNum%binNum;

	if(binSize<1)
	{
		std::cout << "Error for bin size calculation!" << std::endl;
		exit(1);
	}

	int totalVertexOnBin = binSize+1, index = 0;
	float dotValue, leftNorm, rightNorm, meanRotation = 0.0, result, rotationSum;
	Vector3f left, right;
	for (int i = 0; i < binNum; ++i)
	{
		/* would reduce that to binSize if i>=residueNum */
		if(i==residueNum)
			totalVertexOnBin = binSize;

		/* reset the rotationSum */
		rotationSum = 0.0;
		for(int j=0;j<totalVertexOnBin;++j)
		{
			left << array(index*3+3)-array(index*3),
					array(index*3+4)-array(index*3+1),
					array(index*3+5)-array(index*3+2);
			right << array(index*3+6)-array(index*3+3),
					 array(index*3+7)-array(index*3+4),
					 array(index*3+8)-array(index*3+5);
			dotValue = left.dot(right);
			leftNorm = left.norm();
			rightNorm = right.norm();
			if(leftNorm >= 1.0e-8 && rightNorm >=1.0e-8)
			{
				result = dotValue/leftNorm/rightNorm;
				result = min(1.0,(double)result);
				result = max(-1.0,(double)result);
				rotationSum += acos(result);
			}
			++index;
		}

		histogram[i] = rotationSum;
	}
	assert(index==segmentNum);
}
