#include "RANSACModel.h"
#include <vector>



class RANSAC 
{
	typedef RANSACModel::Ptr RANSACModelPtr;

public:
	typedef boost::shared_ptr<RANSAC> Ptr;
	typedef boost::shared_ptr<const RANSAC> ConstPtr;

	int iterations_;
	RANSACModelPtr sac_model_;
	std::vector<int> model_;
	Eigen::VectorXf model_coefficients_;
	std::vector<int> inliers_;

	/** \brief RANSAC (RAndom SAmple Consensus) main constructor
	* \param[in] model a Sample Consensus model
	*/
	RANSAC(const RANSACModelPtr &model)		
		:rng_alg_()
		,rng_(new boost::uniform_01<boost::mt19937>(rng_alg_))
	{
		sac_model_ = model;
		probability_ = 0.99;
		iterations_ = 0;
		threshold_ = std::numeric_limits<double>::max();

		// Maximum number of trials before we give up.
		max_iterations_ = 10000;
	}

	/** \brief RANSAC (RAndom SAmple Consensus) main constructor
	* \param[in] model a Sample Consensus model
	* \param[in] threshold distance to model threshold
	*/
	RANSAC(const RANSACModelPtr &model, double threshold)
		: RANSAC(model, threshold)
	{
		probability_ = 0.99;
		// Maximum number of trials before we give up.
		max_iterations_ = 10000;
	}

	/** \brief Compute the actual model and find the inliers
	* \param[in] debug_verbosity_level enable/disable on-screen debug information and set the verbosity level
	*/
	bool computeModel(int debug_verbosity_level = 0);
	
	bool refineModel(const double sigma = 3.0, const unsigned int max_iterations = 1000);

	inline void	setDistanceThreshold(double threshold) 
	{
		threshold_ = threshold;
	}

	inline double getDistanceThreshold() { return (threshold_); }

	inline void setMaxIterations(int max_iterations) {
		max_iterations_ = max_iterations;
	}

	inline int getMaxIterations() {
		return (max_iterations_);
	}

	inline void setProbability(double probability) { probability_ = probability; }

	inline void getModelCoefficients(Eigen::VectorXf &model_coefficients) {
		model_coefficients = model_coefficients_;
	}

	inline double getProbability() {
		return (probability_);
	}

protected:
	double threshold_;
	int max_iterations_;
	double probability_;
	boost::mt19937 rng_alg_;
	boost::shared_ptr<boost::uniform_01<boost::mt19937> > rng_;

	inline double
		rnd()
	{
		return ((*rng_) ());
	}

};
