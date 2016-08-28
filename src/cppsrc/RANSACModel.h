#ifndef PCL_RANSAC_IMPL_RANSAC_REGISTRATION_SCALE_H_
#define PCL_RANSAC_IMPL_RANSAC_REGISTRATION_SCALE_H_ 

#include <Eigen/Eigen>
#include <Eigen/SVD>
#include <vector>
#include <map>
#include <boost/random.hpp>
#include <boost/exception_ptr.hpp>
#include <set>


using namespace std;

typedef struct
{
	float x;
	float y;
	float z;
}ptXYZ;

typedef struct
{
	float x;
	float y;
	float z;
	unsigned char R;
	unsigned char G;
	unsigned char B;
}ptXYZRGB;

class RANSACModel
{

public:

	typedef boost::shared_ptr<RANSACModel> Ptr;
	typedef boost::shared_ptr<const std::vector<ptXYZ> > PtsConstPtr;
	PtsConstPtr pts; 
	Eigen::Matrix3f m_rot;
	Eigen::Vector3f m_tran;
	double m_scale;

	RANSACModel(const PtsConstPtr &pts,
		bool random = true);

	~RANSACModel(){};

	/** \brief Provide a pointer to the input dataset
	* \param[in] pts the const boost shared pointer to a PointCloud message
	*/
	 void	setInputCloud(const PtsConstPtr &pts);
	
	/** \brief Set the input point pts target.
	* \param[in] target the input poininput_t pts target
	*/
	 void	setInputTarget(const PtsConstPtr &target);

	 inline boost::shared_ptr <std::vector<int> >
		 getIndices() const { return (indices_); }
	 
	 inline unsigned int
		 getSampleSize() const
	 {
		 return sample_size_;
	 }

	 inline void
		 getRandomSamples(const boost::shared_ptr <std::vector<int> > &indices,
		 size_t nr_samples,
		 std::vector<int> &indices_subset)
	 {
		 indices_subset.clear();
		 while (indices_subset.size() < nr_samples)
		 {
			 int a = static_cast<int> (rnd() % (indices->size() - 1));
		 //indices_subset.insert ((*indices)[(int) (indices->size () * (rand () / (RAND_MAX + 1.0)))]);
			 indices_subset.push_back((*indices)[a]);
			 for (int ii = 0; ii < indices_subset.size(); ii++)
			 {
				 int num = 0;
				 if (a == indices_subset[ii])
					 num++;
				 if (num > 1)
					 indices_subset.erase(indices_subset.end());
			 }
		 }
	 }

	 void getSamples(int &iterations, std::vector<int> &samples);
	

	 inline void
		 drawIndexSample(std::vector<int> &sample)
	 {
		 size_t sample_size = sample.size();
		 size_t index_size = shuffled_indices_.size();
		 for (unsigned int i = 0; i < sample_size; ++i)
			 std::swap(shuffled_indices_[i], shuffled_indices_[i + (rnd() % (index_size - i))]);
		 std::copy(shuffled_indices_.begin(), shuffled_indices_.begin() + sample_size, sample.begin());
	 }

	/** \brief Set the input point pts target.
	* \param[in] target the input point pts target
	* \param[in] indices_tgt a vector of point indices to be used from \a target
	*/
	// void	setInputTarget(const vector<ptXYZ> &target, vector<int> &indices_tgt);

	/** \brief Compute a 4x4 rigid transformation matrix from the samples given
	* \param[in] samples the indices found as good candidates for creating a valid model
	* \param[out] model_coefficients the resultant model coefficients
	*/
	bool computeModelCoefficients(const std::vector<int> &samples, Eigen::VectorXf &model_coefficients);


	/** \brief Compute all distances from the transformed points to their correspondences
	* \param[in] model_coefficients the 4x4 transformation matrix
	* \param[out] distances the resultant estimated distances
	*/
	void getDistancesToModel(const Eigen::VectorXf &model_coefficients,
		std::vector<double> &distances);

	/** \brief Select all the points which respect the given model coefficients as inliers.
	* \param[in] model_coefficients the 4x4 transformation matrix
	* \param[in] threshold a maximum admissible distance threshold for determining the inliers from the outliers
	* \param[out] inliers the resultant model inliers
	*/
	void
		selectWithinDistance(const Eigen::VectorXf &model_coefficients,
		const double threshold,
		std::vector<int> &inliers);

	/** \brief Count all the points which respect the given model coefficients as inliers.
	*
	* \param[in] model_coefficients the coefficients of a model that we need to compute distances to
	* \param[in] threshold maximum admissible distance threshold for determining the inliers from the outliers
	* \return the resultant number of inliers
	*/
	virtual int countWithinDistance(const Eigen::VectorXf &model_coefficients,
		const double threshold);

	/** \brief Recompute the 4x4 transformation using the given inlier set
	* \param[in] inliers the data inliers found as supporting the model
	* \param[in] model_coefficients the initial guess for the optimization
	* \param[out] optimized_coefficients the resultant recomputed transformation
	*/
	void optimizeModelCoefficients(const std::vector<int> &inliers,
		const Eigen::VectorXf &model_coefficients, Eigen::VectorXf &optimized_coefficients);

	double computeVariance(const std::vector<double> &error_sqr_dists)
	{
		std::vector<double> dists(error_sqr_dists);
		const size_t medIdx = dists.size() >> 1;
		std::nth_element(dists.begin(), dists.begin() + medIdx, dists.end());
		double median_error_sqr = dists[medIdx];
		return (2.1981 * median_error_sqr);
	}

	double	computeVariance()
	{
		if (error_sqr_dists_.empty())
		{
			return (std::numeric_limits<double>::quiet_NaN());
		}
		return (computeVariance(error_sqr_dists_));
	}


protected:

	unsigned int sample_size_;
	unsigned int model_size_;
	std::vector<double> error_sqr_dists_;;
	boost::shared_ptr<boost::variate_generator< boost::mt19937&, boost::uniform_int<> > > rng_gen_;
	boost::shared_ptr<boost::uniform_int<> > rng_dist_;
	boost::mt19937 rng_alg_;

	bool isSampleGood(const std::vector<int> &samples) const;

	void
		drawIndexSampleRadius(std::vector<int> &sample)
	{
		size_t sample_size = sample.size();
		size_t index_size = shuffled_indices_.size();

		std::swap(shuffled_indices_[0], shuffled_indices_[0 + (rnd() % (index_size - 0))]);
		//const PointT& pt0 = (*input_)[shuffled_indices_[0]];

		std::vector<int> indices;
		std::vector<float> sqr_dists;

		// If indices have been set when the search object was constructed,
		// radiusSearch() expects an index into the indices vector as its
		// first parameter. This can't be determined efficiently, so we use
		// the point instead of the index.
		// Returned indices are converted automatically.
		radiusSearch(input_->at(shuffled_indices_[0]),
			samples_radius_, indices, sqr_dists);

		if (indices.size() < sample_size - 1)
		{
			// radius search failed, make an invalid model
			for (unsigned int i = 1; i < sample_size; ++i)
				shuffled_indices_[i] = shuffled_indices_[0];
		}
		else
		{
			for (unsigned int i = 0; i < sample_size - 1; ++i)
				std::swap(indices[i], indices[i + (rnd() % (indices.size() - i))]);
			for (unsigned int i = 1; i < sample_size; ++i)
				shuffled_indices_[i] = indices[i - 1];
		}

		std::copy(shuffled_indices_.begin(), shuffled_indices_.begin() + sample_size, sample.begin());
	}



	boost::shared_ptr <std::vector<int> > indices_;

	 void	computeSampleDistanceThreshold(const PtsConstPtr &pts);

	 void computeSampleDistanceThreshold(const PtsConstPtr &pts,
		const std::vector<int> &indices);

	void estimateTransformation(const std::vector <ptXYZ> &pts_src, const std::vector<int> &indices_src,
		const std::vector <ptXYZ> &pts_tgt, const std::vector<int> &indices_tgt, Eigen::VectorXf &transform);

	void computeOriginalIndexMapping();

	inline int	rnd()
	{
		return ((*rng_gen_) ());
	}

	int	radiusSearch(const ptXYZ &point, double radius, std::vector<int>& k_indices,
		std::vector<float>& k_sqr_distances, unsigned int max_nn = 0);

	PtsConstPtr target_;
	PtsConstPtr input_;
	std::vector<int> shuffled_indices_;
	boost::shared_ptr <std::vector<int> > indices_tgt_;
	map<int, int> correspondences_;
	double sample_dist_thresh_;
	static const unsigned int max_sample_checks_ = 1000;
	double samples_radius_;
};






#endif //PCL_RANSAC_IMPL_RANSAC_REGISTRATION_SCALE_H
