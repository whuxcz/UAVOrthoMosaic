#include <Eigen/Eigen>
#include <Eigen/SVD>
#include "RANSACModel.h"
#include <vector>
#include <flann/flann.hpp>
#include <cmath>
#include <iostream>
using namespace std;


int computeCentroid(const std::vector < ptXYZ > &pts, Eigen::Vector4f &centroid)
{
	if (pts.empty())
		return (0);

	// Initialize to 0
	centroid.setZero();
	for (size_t i = 0; i < pts.size(); ++i)
	{
		centroid[0] += pts[i].x;
		centroid[1] += pts[i].y;
		centroid[2] += pts[i].z;
	}
	centroid /= static_cast<float> (pts.size());
	centroid[3] = 1;

	return (static_cast<unsigned int> (pts.size()));
}

void demeanPts(const std::vector <ptXYZ> &pts_in, const Eigen::Vector4f &centroid,
	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> &pts_out)
{
	size_t npts = pts_in.size();

	pts_out = Eigen::Matrix<float, 4, Eigen::Dynamic>::Zero(4, npts);        // keep the data aligned

	for (size_t i = 0; i < npts; ++i)
	{
		pts_out(0, i) = pts_in[i].x - centroid[0];
		pts_out(1, i) = pts_in[i].y - centroid[1];
		pts_out(2, i) = pts_in[i].z - centroid[2];
	}
}

 void computeRoots2(const float& b, const float& c, Eigen::Vector3f &roots)
{
	roots(0) = float(0);
	float d = float(b * b - 4.0 * c);
	if (d < 0.0)  // no real roots ! THIS SHOULD NOT HAPPEN!
		d = 0.0;

	float sd = ::sqrt(d);

	roots(2) = 0.5f * (b + sd);
	roots(1) = 0.5f * (b - sd);
}

 unsigned int computeMeanAndCovarianceMatrix(const std::vector<ptXYZ> &pts,
	const std::vector<int> &indices,
	Eigen::Matrix<float, 3, 3> &covariance_matrix,
	Eigen::Matrix<float, 4, 1> &centroid)
{
	// create the buffer on the stack which is much faster than using pts[indices[i]] and centroid as a buffer
	Eigen::Matrix<float, 1, 9, Eigen::RowMajor> accu = Eigen::Matrix<float, 1, 9, Eigen::RowMajor>::Zero();
	size_t point_count;
	point_count = indices.size();
	for (std::vector<int>::const_iterator iIt = indices.begin(); iIt != indices.end(); ++iIt)
	{
		//const ptXYZ& point = pts[*iIt];
		accu[0] += pts[*iIt].x * pts[*iIt].x;
		accu[1] += pts[*iIt].x * pts[*iIt].y;
		accu[2] += pts[*iIt].x * pts[*iIt].z;
		accu[3] += pts[*iIt].y * pts[*iIt].y;
		accu[4] += pts[*iIt].y * pts[*iIt].z;
		accu[5] += pts[*iIt].z * pts[*iIt].z;
		accu[6] += pts[*iIt].x;
		accu[7] += pts[*iIt].y;
		accu[8] += pts[*iIt].z;
	}

	accu /= static_cast<float> (point_count);
	//Eigen::Vector3f vec = accu.tail<3> ();
	//centroid.head<3> () = vec;//= accu.tail<3> ();
	//centroid.head<3> () = accu.tail<3> ();    -- does not compile with Clang 3.0
	centroid[0] = accu[6]; centroid[1] = accu[7]; centroid[2] = accu[8];
	centroid[3] = 1;
	covariance_matrix.coeffRef(0) = accu[0] - accu[6] * accu[6];
	covariance_matrix.coeffRef(1) = accu[1] - accu[6] * accu[7];
	covariance_matrix.coeffRef(2) = accu[2] - accu[6] * accu[8];
	covariance_matrix.coeffRef(4) = accu[3] - accu[7] * accu[7];
	covariance_matrix.coeffRef(5) = accu[4] - accu[7] * accu[8];
	covariance_matrix.coeffRef(8) = accu[5] - accu[8] * accu[8];
	covariance_matrix.coeffRef(3) = covariance_matrix.coeff(1);
	covariance_matrix.coeffRef(6) = covariance_matrix.coeff(2);
	covariance_matrix.coeffRef(7) = covariance_matrix.coeff(5);

	return (static_cast<unsigned int> (point_count));
}

 unsigned int computeMeanAndCovarianceMatrix(const std::vector<ptXYZ> &pts,
Eigen::Matrix<float, 3, 3> &covariance_matrix,
Eigen::Matrix<float, 4, 1> &centroid)
{
	// create the buffer on the stack which is much faster than using pts[indices[i]] and centroid as a buffer
	Eigen::Matrix<float, 1, 9, Eigen::RowMajor> accu = Eigen::Matrix<float, 1, 9, Eigen::RowMajor>::Zero();
	size_t point_count;
	point_count = pts.size();
	// For each point in the pts
	for (size_t i = 0; i < point_count; ++i)
	{
		accu[0] += pts[i].x * pts[i].x;
		accu[1] += pts[i].x * pts[i].y;
		accu[2] += pts[i].x * pts[i].z;
		accu[3] += pts[i].y * pts[i].y; // 4
		accu[4] += pts[i].y * pts[i].z; // 5
		accu[5] += pts[i].z * pts[i].z; // 8
		accu[6] += pts[i].x;
		accu[7] += pts[i].y;
		accu[8] += pts[i].z;
	}
	
	accu /= static_cast<float> (point_count);
	if (point_count != 0)
	{
		//centroid.head<3> () = accu.tail<3> ();    -- does not compile with Clang 3.0
		centroid[0] = accu[6]; centroid[1] = accu[7]; centroid[2] = accu[8];
		centroid[3] = 1;
		covariance_matrix.coeffRef(0) = accu[0] - accu[6] * accu[6];
		covariance_matrix.coeffRef(1) = accu[1] - accu[6] * accu[7];
		covariance_matrix.coeffRef(2) = accu[2] - accu[6] * accu[8];
		covariance_matrix.coeffRef(4) = accu[3] - accu[7] * accu[7];
		covariance_matrix.coeffRef(5) = accu[4] - accu[7] * accu[8];
		covariance_matrix.coeffRef(8) = accu[5] - accu[8] * accu[8];
		covariance_matrix.coeffRef(3) = covariance_matrix.coeff(1);
		covariance_matrix.coeffRef(6) = covariance_matrix.coeff(2);
		covariance_matrix.coeffRef(7) = covariance_matrix.coeff(5);
	}
	return (static_cast<unsigned int> (point_count));
}


void computeRoots(const Eigen::Matrix3f &m, Eigen::Vector3f &roots)
{

	// The characteristic equation is x^3 - c2*x^2 + c1*x - c0 = 0.  The
	// eigenvalues are the roots to this equation, all guaranteed to be
	// real-valued, because the matrix is symmetric.
	float c0 = m(0, 0) * m(1, 1) * m(2, 2)
		+ float(2) * m(0, 1) * m(0, 2) * m(1, 2)
		- m(0, 0) * m(1, 2) * m(1, 2)
		- m(1, 1) * m(0, 2) * m(0, 2)
		- m(2, 2) * m(0, 1) * m(0, 1);
	float c1 = m(0, 0) * m(1, 1) -
		m(0, 1) * m(0, 1) +
		m(0, 0) * m(2, 2) -
		m(0, 2) * m(0, 2) +
		m(1, 1) * m(2, 2) -
		m(1, 2) * m(1, 2);
	float c2 = m(0, 0) + m(1, 1) + m(2, 2);

	if (fabs(c0) < Eigen::NumTraits < float > ::epsilon())  // one root is 0 -> quadratic equation
		computeRoots2(c2, c1, roots);
	else
	{
		const float s_inv3 = float(1.0 / 3.0);
		const float s_sqrt3 = sqrt(float(3.0));
		// Construct the parameters used in classifying the roots of the equation
		// and in solving the equation for the roots in closed form.
		float c2_over_3 = c2 * s_inv3;
		float a_over_3 = (c1 - c2 * c2_over_3) * s_inv3;
		if (a_over_3 > float(0))
			a_over_3 = float(0);

		float half_b = float(0.5) * (c0 + c2_over_3 * (float(2) * c2_over_3 * c2_over_3 - c1));

		float q = half_b * half_b + a_over_3 * a_over_3 * a_over_3;
		if (q > float(0))
			q = float(0);

		// Compute the eigenvalues by solving for the roots of the polynomial.
		float rho = sqrt(-a_over_3);
		float theta = atan2(sqrt(-q), half_b) * s_inv3;
		float cos_theta = cos(theta);
		float sin_theta = sin(theta);
		roots(0) = c2_over_3 + float(2) * rho * cos_theta;
		roots(1) = c2_over_3 - rho * (cos_theta + s_sqrt3 * sin_theta);
		roots(2) = c2_over_3 - rho * (cos_theta - s_sqrt3 * sin_theta);

		// Sort in increasing order.
		if (roots(0) >= roots(1))
			swap(roots(0), roots(1));
		if (roots(1) >= roots(2))
		{
			swap(roots(1), roots(2));
			if (roots(0) >= roots(1))
				swap(roots(0), roots(1));
		}

		if (roots(0) <= 0)  // eigenval for symetric positive semi-definite matrix can not be negative! Set it to 0
			computeRoots2(c2, c1, roots);
	}
}

RANSACModel::RANSACModel(const PtsConstPtr &pts, bool random)
	:input_()
	, indices_()
	, target_()
	, indices_tgt_()
	, correspondences_()
	, sample_dist_thresh_(0)
	, error_sqr_dists_()
	, rng_alg_()
	, rng_dist_(new boost::uniform_int<>(0, std::numeric_limits<int>::max()))
	, rng_gen_()
	, samples_radius_(0.)
{
	// Create a random number generator object
	if (random)
		rng_alg_.seed(static_cast<unsigned> (std::time(0)));
	else
		rng_alg_.seed(12345u);

	rng_gen_.reset(new boost::variate_generator<boost::mt19937&, boost::uniform_int<> >(rng_alg_, *rng_dist_));

	// Call our own setInputCloud
	setInputCloud(pts);
	sample_size_ = 3;
	model_size_ = 16;
	m_rot.setIdentity();
	m_tran.setZero();
	m_scale = 1;
}


 void RANSACModel::setInputCloud(const PtsConstPtr &pts)
{
	input_ = pts;
	if (!indices_)
		indices_.reset(new std::vector<int>());
	if (indices_->empty())
	{
		// Prepare a set of indices to be used (entire pts)
		indices_->resize(pts->size());
		for (size_t i = 0; i < pts->size(); ++i)
			(*indices_)[i] = static_cast<int> (i);
	}
	shuffled_indices_ = *indices_;
	computeOriginalIndexMapping();
	computeSampleDistanceThreshold(pts);
}

void	RANSACModel::setInputTarget(const PtsConstPtr &target)
{
	target_ = target;
	indices_tgt_.reset(new std::vector<int>);
	// Cache the size and fill the target indices
	int target_size = static_cast<int> (target->size());
	indices_tgt_->resize(target_size);

	for (int i = 0; i < target_size; ++i)
		(*indices_tgt_)[i] = i;
	computeOriginalIndexMapping();
}


void RANSACModel::getSamples(int &iterations, std::vector<int> &samples)
{
	// We're assuming that indices_ have already been set in the constructor
	if (indices_->size() < getSampleSize())
	{
		// one of these will make it stop :)
		samples.clear();
		iterations = INT_MAX - 1;
		return;
	}

	// Get a second point which is different than the first
	samples.resize(getSampleSize());
	for (unsigned int iter = 0; iter < max_sample_checks_; ++iter)
	{
		// Choose the random indices
		if (samples_radius_ < std::numeric_limits<double>::epsilon())
			RANSACModel::drawIndexSample(samples);
		else
			RANSACModel::drawIndexSampleRadius(samples);

		// If it's a good sample, stop here
		if (isSampleGood(samples))
		{
			return;
		}
	}
	samples.clear();
}


bool RANSACModel::isSampleGood(const std::vector<int> &samples) const
{
	ptXYZ p10;
	p10.x = input_->at(samples[1]).x - input_->at(samples[0]).x;
	p10.y = input_->at(samples[1]).y - input_->at(samples[0]).y;
	p10.z = input_->at(samples[1]).z - input_->at(samples[0]).z;
	ptXYZ p20; 
	p20.x = input_->at(samples[2]).x - input_->at(samples[0]).x;
	p20.y = input_->at(samples[2]).y - input_->at(samples[0]).y;
	p20.z = input_->at(samples[2]).z - input_->at(samples[0]).z;
	ptXYZ p21;
	p21.x= input_->at(samples[2]).x - input_->at(samples[1]).x;
	p21.y= input_->at(samples[2]).y - input_->at(samples[1]).y;
	p21.z= input_->at(samples[2]).z - input_->at(samples[1]).z;

	return ((p10.x * p10.x + p10.y * p10.y + p10.z * p10.z) > sample_dist_thresh_ &&
		(p20.x * p20.x + p20.y * p20.y + p20.z * p20.z) > sample_dist_thresh_ &&
		(p21.x * p21.x + p21.y * p21.y + p21.z * p21.z) > sample_dist_thresh_);
}

bool RANSACModel::computeModelCoefficients(const std::vector<int> &samples, Eigen::VectorXf &model_coefficients)
{
	if (!target_)
		return (false);

	// Need 3 samples
	if (samples.size() != 3)
		return (false);

	std::vector<int> indices_tgt(3);
	for (int i = 0; i < 3; ++i)
		indices_tgt[i] = correspondences_[samples[i]];

	estimateTransformation(*input_, samples, *target_, indices_tgt, model_coefficients);
	return (true);
}


void RANSACModel::getDistancesToModel(const Eigen::VectorXf &model_coefficients,
	std::vector<double> &distances)
{
	if (indices_->size() != indices_tgt_->size())
	{
		distances.clear();
		return;
	}
	if (!target_)
	{
		return;
	}
	// Check if the model is valid given the user constraints
	if (model_coefficients.size() != model_size_)
	{
		distances.clear();
		return;
	}
	distances.resize(indices_->size());

	// Get the 4x4 transformation
	Eigen::Matrix4f transform;
	transform.row(0).matrix() = model_coefficients.segment<4>(0);
	transform.row(1).matrix() = model_coefficients.segment<4>(4);
	transform.row(2).matrix() = model_coefficients.segment<4>(8);
	transform.row(3).matrix() = model_coefficients.segment<4>(12);

	for (size_t i = 0; i < indices_->size(); ++i)
	{
		Eigen::Vector4f pt_src(input_->at((*indices_)[i]).x,
			input_->at((*indices_)[i]).y,
			input_->at((*indices_)[i]).z, 1);
		Eigen::Vector4f pt_tgt(target_->at((*indices_tgt_)[i]).x,
			target_->at((*indices_tgt_)[i]).y,
			target_->at((*indices_tgt_)[i]).z, 1);

		Eigen::Vector4f p_tr(transform * pt_src);
		// Calculate the distance from the transformed point to its correspondence
		// need to compute the real norm here to keep MSAC and friends general
		distances[i] = (p_tr - pt_tgt).norm();
	}
}

void RANSACModel::selectWithinDistance(const Eigen::VectorXf &model_coefficients, const double threshold, std::vector<int> &inliers)
{
	if (indices_->size() != indices_tgt_->size())
	{
		inliers.clear();
		return;
	}
	if (!target_)
	{
		return;
	}

	double thresh = threshold * threshold;

	// Check if the model is valid given the user constraints
	if (model_coefficients.size() != model_size_)
	{
		inliers.clear();
		return;
	}

	int nr_p = 0;
	inliers.resize(indices_->size());
	error_sqr_dists_.resize(indices_->size());

	Eigen::Matrix4f transform;
	transform.row(0).matrix() = model_coefficients.segment<4>(0);
	transform.row(1).matrix() = model_coefficients.segment<4>(4);
	transform.row(2).matrix() = model_coefficients.segment<4>(8);
	transform.row(3).matrix() = model_coefficients.segment<4>(12);

	for (size_t i = 0; i < indices_->size(); ++i)
	{
		Eigen::Vector4f pt_src(input_->at((*indices_)[i]).x,
			input_->at((*indices_)[i]).y,
			input_->at((*indices_)[i]).z, 1);
		Eigen::Vector4f pt_tgt(target_->at((*indices_tgt_)[i]).x,
			target_->at((*indices_tgt_)[i]).y,
			target_->at((*indices_tgt_)[i]).z, 1);

		Eigen::Vector4f p_tr(transform * pt_src);

		float distance = (p_tr - pt_tgt).squaredNorm();
		// Calculate the distance from the transformed point to its correspondence
		if (distance < thresh)
		{
			inliers[nr_p] = (*indices_)[i];
			error_sqr_dists_[nr_p] = static_cast<double> (distance);
			++nr_p;
		}
	}
	inliers.resize(nr_p);
	error_sqr_dists_.resize(nr_p);
}

int RANSACModel::countWithinDistance(const Eigen::VectorXf &model_coefficients, const double threshold)
{
	if (indices_->size() != indices_tgt_->size())
	{
		return (0);
	}
	if (!target_)
	{
		return (0);
	}

	double thresh = threshold * threshold;

	// Check if the model is valid given the user constraints
	if (model_coefficients.size() != model_size_)
		return (0);

	Eigen::Matrix4f transform;
	transform.row(0).matrix() = model_coefficients.segment<4>(0);
	transform.row(1).matrix() = model_coefficients.segment<4>(4);
	transform.row(2).matrix() = model_coefficients.segment<4>(8);
	transform.row(3).matrix() = model_coefficients.segment<4>(12);

	int nr_p = 0;
	for (size_t i = 0; i < indices_->size(); ++i)
	{
		Eigen::Vector4f pt_src(input_->at((*indices_)[i]).x,
			input_->at((*indices_)[i]).y,
			input_->at((*indices_)[i]).z, 1);
		Eigen::Vector4f pt_tgt(target_->at((*indices_tgt_)[i]).x,
			target_->at((*indices_tgt_)[i]).y,
			target_->at((*indices_tgt_)[i]).z, 1);

		Eigen::Vector4f p_tr(transform * pt_src);
		// Calculate the distance from the transformed point to its correspondence
		if ((p_tr - pt_tgt).squaredNorm() < thresh)
			nr_p++;
	}
	return (nr_p);
}


void RANSACModel::optimizeModelCoefficients(const std::vector<int> &inliers, const Eigen::VectorXf &model_coefficients, Eigen::VectorXf &optimized_coefficients)
{
	if (indices_->size() != indices_tgt_->size())
	{
		optimized_coefficients = model_coefficients;
		return;
	}

	// Check if the model is valid given the user constraints
	if (model_coefficients.size() != model_size_ || !target_)
	{
		optimized_coefficients = model_coefficients;
		return;
	}

	std::vector<int> indices_src(inliers.size());
	std::vector<int> indices_tgt(inliers.size());
	for (size_t i = 0; i < inliers.size(); ++i)
	{
		indices_src[i] = inliers[i];
		indices_tgt[i] = correspondences_[indices_src[i]];
	}

	estimateTransformation(*input_, indices_src, *target_, indices_tgt, optimized_coefficients);
}

 void	RANSACModel::computeSampleDistanceThreshold(const PtsConstPtr &pts)
{
	// Compute the principal directions via PCA
	Eigen::Vector4f xyz_centroid;
	Eigen::Matrix3f covariance_matrix = Eigen::Matrix3f::Zero();

	computeMeanAndCovarianceMatrix(*pts, covariance_matrix, xyz_centroid);

	// Check if the covariance matrix is finite or not.
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			if (!std::isfinite(covariance_matrix.coeffRef(i, j)))
				exit(-1);
	Eigen::Vector3f eigen_values;

	float scale = covariance_matrix.cwiseAbs().maxCoeff();
	if (scale <= numeric_limits < float > ::min())
		scale = float(1.0);

	Eigen::Matrix3f scaledMat = covariance_matrix / scale;
	computeRoots(scaledMat, eigen_values);
	eigen_values *= scale;

	// Compute the distance threshold for sample selection
	sample_dist_thresh_ = eigen_values.array().sqrt().sum() / 3.0;
	sample_dist_thresh_ *= sample_dist_thresh_;
}

 void RANSACModel::computeSampleDistanceThreshold(const PtsConstPtr &pts,
	const std::vector<int> &indices)
{
	// Compute the principal directions via PCA
	Eigen::Vector4f xyz_centroid;
	Eigen::Matrix3f covariance_matrix;
	computeMeanAndCovarianceMatrix(*pts, indices, covariance_matrix, xyz_centroid);

	// Check if the covariance matrix is finite or not.
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			if (!std::isfinite(covariance_matrix.coeffRef(i, j)))
				exit(-1);

	Eigen::Vector3f eigen_values;

	float scale = covariance_matrix.cwiseAbs().maxCoeff();

	if (scale <= numeric_limits < float > ::min())
		scale = float(1.0);

	Eigen::Matrix3f scaledMat = covariance_matrix / scale;
	computeRoots(scaledMat, eigen_values);
	eigen_values *= scale;

	// Compute the distance threshold for sample selection
	sample_dist_thresh_ = eigen_values.array().sqrt().sum() / 3.0;
	sample_dist_thresh_ *= sample_dist_thresh_;
}

void RANSACModel::estimateTransformation(const std::vector <ptXYZ> &pts_src, const std::vector<int> &indices_src,
	const std::vector <ptXYZ> &pts_tgt, const std::vector<int> &indices_tgt, Eigen::VectorXf &transform)
{
	transform.resize(16);

	std::vector<ptXYZ> src;
	std::vector<ptXYZ> tgt;

	for (size_t i = 0; i < indices_src.size(); i++)
	{
		src.push_back(pts_src[indices_src[i]]);
		tgt.push_back(pts_tgt[indices_tgt[i]]);
	}

	Eigen::Matrix<float, 4, 1> centroid_src;
	Eigen::Matrix<float, 4, 1> centroid_tgt;
	computeCentroid(src, centroid_src);
	computeCentroid(tgt, centroid_tgt);


	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> pts_src_demean, pts_tgt_demean;
	demeanPts(src, centroid_src, pts_src_demean);
	demeanPts(tgt, centroid_tgt, pts_tgt_demean);
	Eigen::Matrix<float, 4, 4> mtransform;
	mtransform.setIdentity();

	Eigen::Matrix<float, 3, 3> H = (pts_src_demean * pts_tgt_demean.transpose()).topLeftCorner(3, 3);

	Eigen::JacobiSVD<Eigen::Matrix<float, 3, 3> > svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix<float, 3, 3> u = svd.matrixU();
	Eigen::Matrix<float, 3, 3> v = svd.matrixV();

	if (u.determinant() * v.determinant() < 0)
	{
		for (int x = 0; x < 3; ++x)
			v(x, 2) *= -1;
	}

	Eigen::Matrix<float, 3, 3> R = v * u.transpose();

	m_rot = R;

	// rotated pts
	Eigen::Matrix<float, 4, 4> R4;
	R4.block(0, 0, 3, 3) = R;
	R4(0, 3) = 0;
	R4(1, 3) = 0;
	R4(2, 3) = 0;
	R4(3, 3) = 1;

	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> src_ = R4 * pts_src_demean;

	float scale1, scale2;
	double sum_ss = 0.0f, sum_tt = 0.0f, sum_tt_ = 0.0f;
	for (unsigned corrIdx = 0; corrIdx < pts_src_demean.cols(); ++corrIdx)
	{
		sum_ss += pts_src_demean(0, corrIdx) * pts_src_demean(0, corrIdx);
		sum_ss += pts_src_demean(1, corrIdx) * pts_src_demean(1, corrIdx);
		sum_ss += pts_src_demean(2, corrIdx) * pts_src_demean(2, corrIdx);

		sum_tt += pts_tgt_demean(0, corrIdx) * pts_tgt_demean(0, corrIdx);
		sum_tt += pts_tgt_demean(1, corrIdx) * pts_tgt_demean(1, corrIdx);
		sum_tt += pts_tgt_demean(2, corrIdx) * pts_tgt_demean(2, corrIdx);

		sum_tt_ += pts_tgt_demean(0, corrIdx) * src_(0, corrIdx);
		sum_tt_ += pts_tgt_demean(1, corrIdx) * src_(1, corrIdx);
		sum_tt_ += pts_tgt_demean(2, corrIdx) * src_(2, corrIdx);
	}

	scale1 = sqrt(sum_tt / sum_ss);
	scale2 = sum_tt_ / sum_ss;
	float scale = scale2;
	mtransform.topLeftCorner(3, 3) = scale * R;
	const Eigen::Matrix<float, 3, 1> Rc(scale * R * centroid_src.head(3));
	mtransform.block(0, 3, 3, 1) = centroid_tgt.head(3) - Rc;
	transform.segment<4>(0) = mtransform.row(0).matrix();
	transform.segment<4>(4) = mtransform.row(1).matrix();
	transform.segment<4>(8) = mtransform.row(2).matrix();
	transform.segment<4>(12) = mtransform.row(3).matrix();
	m_tran = centroid_tgt.head(3) - Rc;
	m_scale = scale;

}



void RANSACModel::computeOriginalIndexMapping()
{
	if (!indices_tgt_ || !indices_ || indices_->empty() || indices_->size() != indices_tgt_->size())
		return;
	for (size_t i = 0; i < indices_->size(); ++i)
		correspondences_[(*indices_)[i]] = (*indices_tgt_)[i];
}


int	RANSACModel::radiusSearch(const ptXYZ &point, double radius, std::vector<int>& k_indices,
	std::vector<float>& k_sqr_distances, unsigned int max_nn)
{
	int num_pts = input_->size();
	float* p = new float[input_->size()*3];
	for (int i = 0; i < input_->size(); i++)
	{
		*(p + i) = input_->at(i).x;
		*(p + i + 1) = input_->at(i).y;
		*(p + i + 2) = input_->at(i).z;
	}
	flann::Matrix<float> features(p, num_pts, 3);
	flann::Index<flann::L2<float> > index(features, flann::KDTreeSingleIndexParams());
	index.buildIndex();
	float pt[3] = { point.x, point.y ,point.z};
	std::vector<std::vector<int> > indices;
	std::vector<std::vector<float> > dists;
	flann::Matrix<float> query(pt, 1, 3);
	int count = index.radiusSearch(query, indices,
		dists, radius, flann::SearchParams(10));
	for (size_t i = 0; i < indices.size(); i++)
	{
		k_indices.push_back(indices[0][i]);
		k_sqr_distances.push_back(dists[0][i]);
	}
	return indices.size();
}
