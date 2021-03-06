/*
 * @brief The class contains some preliminary functions to calculate before the distance matrix starts
 * @author Lieyu Shi
 */


#include "PreComputing.h"


/*
 * @brief calculate the sequence values for a given streamline coordinate
 *
 * @param[in] array The coordinate of a streamline
 * @param[in] size The size of coordinate
 * @param[out] rowSequence The rowSequence vector
 */
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


/*
 * @brief Calculate the average discrete curvatures on the streamline coordinate
 *
 * @param[in] array The coordinate of the streamline
 * @param[in] size The size of coordinate
 * @return The float value of the average rotation
 */
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


/*
 * @brief Get the normalized multivariate w.r.t. segments of the streamline
 *
 * @param[in] array The coordinate of the streamline
 * @param[in] size The size of coordinate
 * @param[out] rowSequence The row sequence to be updated
 */
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


/*
 * @brief Get the fixed sequence for the streamline
 *
 * @param[in] array The coordinate of the streamline
 * @param[in] The size of coordinate
 * @param[out] rowSequence The row sequence to be updated
 */
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
			angle = M_PI;
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


/*
 * @brief Get the unnormalized multivariate w.r.t. segments of the streamline coordinates
 *
 * @param[in] array The coordinate of the streamline
 * @param[in] The size of coordinate
 * @param[out] rowSequence The row sequence to be updated
 */
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


/*
 * @brief Get the unit direction for each streamline
 *
 * @param[in] array The coordinate of a streamline
 * @param[in] pointNum The size of point in the line
 * @param[out] direction The direction vector to be updated
 */
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


/*
 * @brief Get the pairwise each attribute for the streamline
 *
 * @param[in] array The coordinate of a streamline
 * @param[in] size The size of the coordinate
 * @param[out] wiseVec The pairwise intersection angle of each segment
 * @param[out] wiseNorm The pairwise norm vector of each segment
 */
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


/*
 * @brief The difference compared to the former function is that, this will resample on high-curvature points
 *
 * @param[in] array The coordinate of a streamline
 * @param[in] binNum The number of points on the line
 * @param[out] histogram The vector values making up the histogram to be updated
 */
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


/*
 * @brief Get the bin-based histogram for signature, the difference is that we should try to get maximal binNum-1 points
 *
 * @param[in] array The coordinate of the streamline
 * @param[in] binNum The number of bins on the histogram
 * @param[out] histogram The histogram of values as an output of vector
 */
void getSignatureHistSampled(const Eigen::VectorXf& array,
					  	  	 const int& binNum,
							 std::vector<float>& histogram)
{
	/* if empty vector, should allocate memory ahead of time */
	if(histogram.empty())
		histogram = std::vector<float>(binNum);

	/* get how many vertices you'll have */
	const int& segmentNum = array.size()/3-2;

	/* preset a priority_queue to get the sampled points in maximal curvatures */
	priority_queue<CurvatureObject, std::vector<CurvatureObject>, CompareFunc> pQueue;

	std::vector<float> curvatureVec(segmentNum);

	Eigen::Vector3f firstSeg, secondSeg;

	/* discrete curvature */
	float curva;

	int vecIndex = 0;
	for(int i=0;i<segmentNum;++i)
	{
		for(int j=0;j<3;++j)
		{
			firstSeg(j)=array(3*i+3+j)-array(3*i+j);
			secondSeg(j)=array(3*i+6+j)-array(3*i+3+j);
		}

		float firstSegNorm = firstSeg.norm(), secondSegNorm = secondSeg.norm();
		if(firstSegNorm<1.0e-8 || secondSegNorm<1.0e-8)
			curva = 0.0;
		else
		{
			curva = firstSeg.dot(secondSeg)/firstSegNorm/secondSegNorm;

			/* clip curvature into range [-1.0, 1.0] */
			curva = std::min(float(1.0), curva);
			curva = std::max(float(-1.0), curva);

			curva = acos(curva);
		}
		/* store in the vector */
		curvatureVec[vecIndex++]=curva;

		/* push it into the priority queue */
		pQueue.push(CurvatureObject(curva, i));

	}

	/* get the first binNum-1 object */
	CurvatureObject top;

	/* use ordered_set to sort the index */
	std::vector<int> indexVec;

	int indexNum = 0;
	const int& requiredNum = binNum-1;
	while(indexNum<requiredNum && !pQueue.empty())
	{
		top = pQueue.top();
		indexVec.push_back(top.index);
		pQueue.pop();
		++indexNum;
	}

	assert(indexVec.size()==requiredNum);

	/* sort the vec */
	std::sort(indexVec.begin(), indexVec.end());

	/* start sampling to make a curvature histogram */
	float curvatureSum = 0.0;

	/* get accumulative curvature */
	int left = 0, right;
	for(int i=0;i<requiredNum;++i)
	{
		right = indexVec[i];

		/* sum up the curvature of left and right */
		curvatureSum = 0.0;
		for(int j=left;j<=right;++j)
		{
			curvatureSum+=curvatureVec[j];
		}

		histogram[i] = curvatureSum;

		left = right+1;
	}

	/* add last element which is from left to last vertex */
	curvatureSum = 0.0;
	for(int i=left;i<segmentNum;++i)
		curvatureSum+=curvatureVec[i];
	histogram[requiredNum] = curvatureSum;
}


/*
 * @brief Get the linear and angular entropy
 *
 * @param[in] array The coordinate of the streamline
 * @param[in] bundleSize The size of bundles for the attributes generated
 * @param[out] histogram The histogram of values as output
 */
void getLinearAngularEntropy(const Eigen::VectorXf& array,
 	  	 	 	 	 	 	 const int& bundleSize,
							 std::vector<float>& histogram)
{
	/* if empty vector, should allocate memory ahead of time */
	if(histogram.empty())
		histogram = std::vector<float>(2);

	/* get how many vertices you'll have */
	const int& segmentNum = array.size()/3-1;

	const int& curvatureNum = segmentNum-1;

	/* should partition the whole streamlines into bunleSize segments, and compute the entropy */

	std::vector<float> segmentVec(segmentNum), curvatureVec(curvatureNum);

	Eigen::Vector3f firstSeg, secondSeg;

	/* discrete curvature */
	float curva;

	float lengthSum = 0.0, curveSum = 0.0;
	int vecIndex = 0;
	for(int i=0;i<curvatureNum;++i)
	{
		for(int j=0;j<3;++j)
		{
			firstSeg(j)=array(3*i+3+j)-array(3*i+j);
			secondSeg(j)=array(3*i+6+j)-array(3*i+3+j);
		}

		if(firstSeg.norm()<1.0e-6 || secondSeg.norm()<1.0e-6)
		{
			segmentVec[i] = 0.0;
			curvatureVec[i]= 0.0;
			continue;
		}

		float firstSegNorm = firstSeg.norm(), secondSegNorm = secondSeg.norm();
		if(firstSegNorm<1.0e-8 || secondSegNorm<1.0e-8)
			curva = 0.0;
		else
		{
			curva = firstSeg.dot(secondSeg)/firstSegNorm/secondSegNorm;

			/* clip curvature into range [-1.0, 1.0] */
			curva = std::min(float(1.0), curva);
			curva = std::max(float(-1.0), curva);

			curva = acos(curva);
		}

		/* store in the vector */
		curvatureVec[i]=curva;
		curveSum+=curva;

		/* store path */
		segmentVec[i] = firstSeg.norm();
		lengthSum+=segmentVec[i];
	}

	int i = curvatureNum;
	for(int j=0;j<3;++j)
	{
		firstSeg(j)=array(3*i+3+j)-array(3*i+j);
	}
	segmentVec[i] = firstSeg.norm();
	lengthSum+=segmentVec[i];


	/* should deal with exceptional case if lengthSum == 0 or curveSum == 0 */
	if(lengthSum<1.0e-6)
	{
		histogram[0] = 1.0;
	}
	else
	{
		/* get ratio for the vec */
		const int& segmentQuotient = segmentNum/bundleSize;
		const int& segmentResidue = segmentNum%bundleSize;

		/* get the vec for bundleSize */
		std::vector<float> lengthVec(bundleSize);

		float tempLength, linearEntropy = 0.0, prob;
		int left, right;
		for(int k = 0;k<bundleSize-1;++k)
		{
			tempLength = 0.0;
			left = k*segmentQuotient, right = (k+1)*segmentQuotient;
			for(int i = left;i<right;++i)
				tempLength+=segmentVec[i];

			prob = tempLength/lengthSum;

			if(prob>1.0e-6)
				linearEntropy += prob*log2f(prob);
		}

		left = (bundleSize-1)*segmentQuotient, right = segmentNum;
		tempLength = 0.0;
		for(int i=left;i<right;++i)
		{
			tempLength+=segmentVec[i];
		}
		if(prob>1.0e-6)
			prob = tempLength/lengthSum;
		linearEntropy += prob*log2f(prob);

		linearEntropy = -linearEntropy/log2f(float(bundleSize));
		histogram[0] = linearEntropy;
	}

	/* deal with curveSum == 1.0 */
	if(curveSum<1.0e-6)
	{
		histogram[1] = 1.0;
	}
	else
	{
		const int& curvatureQuotient = curvatureNum/bundleSize;
		const int& curvatureResidue = curvatureNum%bundleSize;

		/* get the vec for bundleSize */
		std::vector<float> curveVec(bundleSize);

		float tempCurve, angularEntropy = 0.0, prob;
		int left, right;
		for(int k = 0;k<bundleSize-1;++k)
		{
			tempCurve = 0.0;
			left = k*curvatureQuotient, right = (k+1)*curvatureQuotient;
			for(int i=left;i<right;++i)
				tempCurve+=curvatureVec[i];

			prob = tempCurve/curveSum;
			if(prob>1.0e-6)
				angularEntropy += prob*log2f(prob);
		}

		left = (bundleSize-1)*curvatureQuotient, right = curvatureNum;
		tempCurve = 0.0;
		for(int i=left;i<right;++i)
		{
			tempCurve+=curvatureVec[i];
		}
		prob = tempCurve/curveSum;
		if(prob>1.0e-6)
			angularEntropy += prob*log2f(prob);

		angularEntropy = -angularEntropy/log2f(float(bundleSize));
		histogram[1] = angularEntropy;
	}
}


