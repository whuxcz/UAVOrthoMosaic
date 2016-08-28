#include <ctime>
#include "RANSACModel.h"
#include "RANSAC.h"
#include <boost/random.hpp>
#include <limits.h>
#include <iostream>


bool RANSAC::computeModel(int)
{
	// Warn and exit if no threshold was set
	if (threshold_ == std::numeric_limits<double>::max())
	{
		return (false);
	}

	iterations_ = 0;
	int n_best_inliers_count = -INT_MAX;
	double k = 1.0;

	std::vector<int> selection;
	Eigen::VectorXf model_coefficients;

	double log_probability = log(1.0 - probability_);
	double one_over_indices = 1.0 / static_cast<double> (sac_model_->getIndices()->size());

	int n_inliers_count = 0;
	unsigned skipped_count = 0;
	// supress infinite loops by just allowing 10 x maximum allowed iterations for invalid model parameters!
	const unsigned max_skip = max_iterations_ * 10;
	// Iterate
	while (iterations_ < k && skipped_count < max_skip)
	{
		// Get X samples which satisfy the model criteria
		//sac_model_->getSamples(iterations_, selection);
		sac_model_->getRandomSamples(sac_model_->getIndices(),3, selection);

		if (selection.empty())
		{
			break;
		}
		
		// Search for inliers in the point cloud for the current plane model M
		if (!sac_model_->computeModelCoefficients(selection, model_coefficients))
		{
			//++iterations_;
			++skipped_count;
			continue;
		}
		
		// Select the inliers that are within threshold_ from the model
		//sac_model_->selectWithinDistance (model_coefficients, threshold_, inliers);
		//if (inliers.empty () && k > 1.0)
		//  continue;

		n_inliers_count = sac_model_->countWithinDistance(model_coefficients, threshold_);
		// Better match ?
		if (n_inliers_count > n_best_inliers_count)
		{
			n_best_inliers_count = n_inliers_count;

			// Save the current model/inlier/coefficients selection as being the best so far
			model_ = selection;
			model_coefficients_ = model_coefficients;

			// Compute the k parameter (k=log(z)/log(1-w^n))
			double w = static_cast<double> (n_best_inliers_count)* one_over_indices;
			double p_no_outliers = 1.0 - pow(w, static_cast<double> (selection.size()));
			p_no_outliers = (std::max) (std::numeric_limits<double>::epsilon(), p_no_outliers);       // Avoid division by -Inf
			p_no_outliers = (std::min) (1.0 - std::numeric_limits<double>::epsilon(), p_no_outliers);   // Avoid division by 0.
			k = log_probability / log(p_no_outliers);
		}

		++iterations_;
		if (iterations_ > max_iterations_)
		{
			break;
		}
	}

	if (model_.empty())
	{
		inliers_.clear();
		return (false);
	}

	// Get the set of inliers that correspond to the best model found so far
	sac_model_->selectWithinDistance(model_coefficients_, threshold_, inliers_);
	return (true);
}

bool RANSAC::refineModel(const double sigma, const unsigned int max_iterations)
{
	if (!sac_model_)
	{
		return (false);
	}

	double inlier_distance_threshold_sqr = threshold_ * threshold_,
		error_threshold = threshold_;
	double sigma_sqr = sigma * sigma;
	unsigned int refine_iterations = 0;
	bool inlier_changed = false, oscillating = false;
	std::vector<int> new_inliers, prev_inliers = inliers_;
	std::vector<size_t> inliers_sizes;
	Eigen::VectorXf new_model_coefficients = model_coefficients_;
	do
	{
		// Optimize the model coefficients
		sac_model_->optimizeModelCoefficients(prev_inliers, new_model_coefficients, new_model_coefficients);
		inliers_sizes.push_back(prev_inliers.size());

		// Select the new inliers based on the optimized coefficients and new threshold
		sac_model_->selectWithinDistance(new_model_coefficients, error_threshold, new_inliers);

		if (new_inliers.empty())
		{
			refine_iterations++;
			if (refine_iterations >= max_iterations)
				break;
			continue;
			//return (false);
		}

		// Estimate the variance and the new threshold
		double variance = sac_model_->computeVariance();
		error_threshold = sqrt(std::min(inlier_distance_threshold_sqr, sigma_sqr * variance));

		inlier_changed = false;
		std::swap(prev_inliers, new_inliers);
		// If the number of inliers changed, then we are still optimizing
		if (new_inliers.size() != prev_inliers.size())
		{
			// Check if the number of inliers is oscillating in between two values
			if (inliers_sizes.size() >= 4)
			{
				if (inliers_sizes[inliers_sizes.size() - 1] == inliers_sizes[inliers_sizes.size() - 3] &&
					inliers_sizes[inliers_sizes.size() - 2] == inliers_sizes[inliers_sizes.size() - 4])
				{
					oscillating = true;
					break;
				}
			}
			inlier_changed = true;
			continue;
		}

		// Check the values of the inlier set
		for (size_t i = 0; i < prev_inliers.size(); ++i)
		{
			// If the value of the inliers changed, then we are still optimizing
			if (prev_inliers[i] != new_inliers[i])
			{
				inlier_changed = true;
				break;
			}
		}
	} while (inlier_changed && ++refine_iterations < max_iterations);

	// If the new set of inliers is empty, we didn't do a good job refining
	if (new_inliers.empty())
	{
		return (false);
	}

	if (oscillating)
	{
		return (true);
	}

	// If no inliers have been changed anymore, then the refinement was successful
	if (!inlier_changed)
	{
		std::swap(inliers_, new_inliers);
		model_coefficients_ = new_model_coefficients;
		return (true);
	}
	return (false);
}


