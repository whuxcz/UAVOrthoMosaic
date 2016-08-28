#include <iomanip>
#include <time.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <vector>
#include <algorithm>
#include <Eigen/SVD>
#include <Eigen/Eigen>
#include <boost/make_shared.hpp>
#include "RANSACModel.h"
#include "RANSAC.h"
#include "exif.h"
#include "coorconv.hpp"
#include "keys2a_flann.h"

typedef struct
{
	double lat;
	double log;
	double alt;
	double f_mm;
	double x;
	double y;
}exifinfo;

typedef struct
{
	double x;
	double y;
}proCent;

int ReadImgList(char* list_in, std::vector<std::string>& img_files) {
	FILE* fp;

	if ((fp = fopen(list_in, "r")) == NULL) {
		printf("Error opening file %s for reading.\n", list_in);
		return 1;
	}

	char buf[512], *start, *end;
	while (fgets(buf, 512, fp)) {
		// Remove trailing new-line
		if (buf[strlen(buf) - 1] == '\n') buf[strlen(buf) - 1] = '\0';

		// Find first non-space character
		start = buf;
		while (isspace(*start)) start++;

		// Skip empty lines
		if (strlen(start) == 0) continue;
		// Append file-name to img_files
		std::string filename(std::strtok(start, " "));
		img_files.push_back(filename);
	}
	fclose(fp);

	// Check we found input files
	if (img_files.size() == 0) {
		printf("No input files found in %s.\n", list_in);
		return 1;
	}

	return 0;
}

int ReadExif(std::string &inImg, exifinfo& info)
{
	FILE *fp = fopen(inImg.c_str(), "rb");

	if (!fp)
	{
		printf("open image: %swrong!\n ", inImg.data());
		return -1;
	}

	fseek(fp, 0, SEEK_END);

	unsigned long fsize = ftell(fp);
	rewind(fp);

	unsigned char *buf = new unsigned char[fsize];
	if (fread(buf, 1, fsize, fp) != fsize)
	{
		printf("Can't read file.\n");
		delete[] buf;
		return -2;
	}

	fclose(fp);

	easyexif::EXIFInfo result;

	int code = result.parseFrom(buf, fsize);

	delete[] buf;
	if (code)
	{
		printf("Error parsing EXIF: code %d\n", code);
		return -3;
	}

	info.lat = result.GeoLocation.Latitude;
	info.log = result.GeoLocation.Longitude;
	info.alt = result.GeoLocation.Altitude;
	info.f_mm = result.FocalLength;

	WGSUTMTransform wut;
	int zone = static_cast<int>((info.log + 180) / 6) + 1;
	UTMCoor uc;
	wut.LatLonToUTMXY(wut.DegToRad(info.lat), wut.DegToRad(info.log), zone, uc);
	info.x = uc.x;
	info.y = uc.y;
	return zone;
}


int GetCamInfo(float *cam_tgt, std::vector<exifinfo> vexif, double &deltaX, double &deltaY)
{
	double sum_x = 0;
	double sum_y = 0;
	proCent center;

	float *tgt = cam_tgt;
	for (size_t i = 0; i < vexif.size(); i++)
	{
		exifinfo info = vexif[i];
		*(tgt + i * 3) = info.x;
		*(tgt + i * 3 + 1) = info.y;
		*(tgt + i * 3 + 2) = info.alt;

		sum_x += info.x;
		sum_y += info.y;
	}
	deltaX = sum_x / vexif.size();
	deltaY = sum_y / vexif.size();

	return 1;
}


int TransformPointCloud(const std::vector<ptXYZ> &pts_in, std::vector<ptXYZ> &pts_out, const Eigen::Matrix4f &transform)
{
	
	if (pts_in.size() != pts_out.size())
	{
		pts_out.clear();
		pts_out.resize(pts_in.size());
	}
	for (size_t i = 0; i < pts_in.size(); i++)
	{
		Eigen::Vector3f pt(pts_in[i].x, pts_in[i].y, pts_in[i].z);

		pts_out[i].x = static_cast<float> (transform(0, 0) * pt.coeffRef(0) + transform(0, 1) * pt.coeffRef(1) + transform(0, 2) * pt.coeffRef(2) + transform(0, 3));
		pts_out[i].y = static_cast<float> (transform(1, 0) * pt.coeffRef(0) + transform(1, 1) * pt.coeffRef(1) + transform(1, 2) * pt.coeffRef(2) + transform(1, 3));
		pts_out[i].z = static_cast<float> (transform(2, 0) * pt.coeffRef(0) + transform(2, 1) * pt.coeffRef(1) + transform(2, 2) * pt.coeffRef(2) + transform(2, 3));
	}

	return pts_in.size();
}


int TransformPointCloud(const std::vector<ptXYZRGB> &pts_in, std::vector<ptXYZRGB> &pts_out, const Eigen::Matrix4f &transform)
{
	if (pts_in.size() != pts_out.size())
	{
		pts_out.clear();
		pts_out.resize(pts_in.size());
	}
	for (size_t i = 0; i < pts_in.size(); i++)
	{
		Eigen::Vector3f pt(pts_in[i].x, pts_in[i].y, pts_in[i].z);

		pts_out[i].x = static_cast<float> (transform(0, 0) * pt.coeffRef(0) + transform(0, 1) * pt.coeffRef(1) + transform(0, 2) * pt.coeffRef(2) + transform(0, 3));
		pts_out[i].y = static_cast<float> (transform(1, 0) * pt.coeffRef(0) + transform(1, 1) * pt.coeffRef(1) + transform(1, 2) * pt.coeffRef(2) + transform(1, 3));
		pts_out[i].z = static_cast<float> (transform(2, 0) * pt.coeffRef(0) + transform(2, 1) * pt.coeffRef(1) + transform(2, 2) * pt.coeffRef(2) + transform(2, 3));
		pts_out[i].R = pts_in[i].R;
		pts_out[i].G = pts_in[i].B;
		pts_out[i].B = pts_in[i].B;
	}

	return pts_in.size();
	return 1;
}

int main(int argc, char* argv[])
{
	if (argc != 6 && argc != 7)
	{
		std::cout << "[Usage] genDem img_list bundle_file dem_file EO.txt pts_utm.ply" << std::endl;
		std::cout << "[OR]    genDem img_list bundle_file dem_file EO.txt pts_utm.ply dense_cloud_file" << std::endl;

		return -1;
	}

	bool is_dense = false;
	if (argc == 7)
		is_dense = true;

	char *list_in;
	char *bundle_in;
	list_in = argv[1];
	bundle_in = argv[2];
	std::vector<std::string> img_files;
	if (ReadImgList(list_in, img_files) != 0)	return -1;

	std::vector<exifinfo> vexif;

	double sumX = 0;
	double sumY = 0;
	std::string line;
	int cam_num = img_files.size();
	float *cam_tgt = new float[cam_num * 3];
	float *cam_src = new float[cam_num * 3];
	float *eos = new float[cam_num * 7];
	int zone;
	for (size_t i = 0; i <img_files.size(); i++)
	{
		if (img_files[i] != "")
		{
			exifinfo info;
			zone = ReadExif(img_files[i], info);
			if (zone)
				vexif.push_back(info);
			else return -1;
		}
	}

	ifstream iBundle;
	iBundle.open(bundle_in);
	int pts_num;

	if (!iBundle)
	{
		std::cout << "Wrong bundle.out file!" << std::endl;
		return -2;
	}
	std::getline(iBundle, line);
	if (line != "# Bundle file v0.3")
	{
		std::cout << "Wrong format of bundle.out!" << std::endl;
		return -3;
	}
	iBundle >> cam_num >> pts_num;
	std::getline(iBundle, line);

	for (size_t i = 0; i < cam_num; i++)
	{
		float f, k1, k2;
		iBundle >> f >> k1 >> k2;
		float r1, r2, r3;
		float r4, r5, r6;
		float r7, r8, r9;
		iBundle >> r1 >> r2 >> r3;
		iBundle >> r4 >> r5 >> r6;
		iBundle >> r7 >> r8 >> r9;
		float t1, t2, t3;
		iBundle >> t1 >> t2 >> t3;
		Eigen::Matrix3f rot;
		Eigen::Vector3f t;
		rot << r1, r2, r3,
			r4, r5, r6,
			r7, r8, r9;
		t << t1, t2, t3;

		double phi = std::atan2(r8, r9);
		double omega = std::asin(r7);
		double kappa = std::atan2(r4, r1);

		rot = rot * -1;
		Eigen::Matrix3f tran_rot = rot.transpose();
		Eigen::Vector3f cam = tran_rot * t;

		*(cam_src + 3 * i) = cam(0);
		*(cam_src + 3 * i + 1) = cam(1);
		*(cam_src + 3 * i + 2) = cam(2);

		*(eos + 7 * i) = phi;
		*(eos + 7 * i + 1) = omega;
		*(eos + 7 * i + 2) = kappa;
		*(eos + 7 * i + 3) = t1;
		*(eos + 7 * i + 4) = t2;
		*(eos + 7 * i + 5) = t3;
		*(eos + 7 * i + 6) = f;
	}
	float minx = std::numeric_limits<float>::max();
	float miny = std::numeric_limits<float>::max();
	float maxx = std::numeric_limits<float>::min();
	float maxy = std::numeric_limits<float>::min();
	std::vector<ptXYZ> pts_src;
	std::vector<ptXYZ> pts_tgt;

	double deltaX = 0;
	double deltaY = 0;
	GetCamInfo(cam_tgt, vexif, deltaX, deltaY);
	for (size_t i = 0; i < cam_num; i++)
	{
		*(cam_tgt + 3 * i) -= deltaX;
		*(cam_tgt + 3 * i + 1) -= deltaY;
	}
	//get transform from model space to UTM
	boost::shared_ptr<std::vector<ptXYZ> > vCamsrc = boost::make_shared<std::vector<ptXYZ> >();
	boost::shared_ptr<std::vector<ptXYZ> > vCamtgt = boost::make_shared<std::vector<ptXYZ> >();
	//std::vector<ptXYZ> vCam_src_;
	std::vector<ptXYZ> vCam_tgt_;

	for (size_t i = 0; i < cam_num; i++)
	{
		ptXYZ pt1, pt2;
		pt1.x = *(cam_src + 3 * i);
		pt1.y = *(cam_src + 3 * i + 1);
		pt1.z = *(cam_src + 3 * i + 2);

		pt2.x = *(cam_tgt + 3 * i);
		pt2.y = *(cam_tgt + 3 * i + 1);
		pt2.z = *(cam_tgt + 3 * i + 2);
		vCam_tgt_.push_back(pt2);
		//vCam_src_.push_back(pt1);
		vCamsrc->push_back(pt1);
		vCamtgt->push_back(pt2);
	}
	

	RANSACModel::Ptr ransac_model(new RANSACModel(vCamsrc));
	ransac_model->setInputTarget(vCamtgt);
	RANSAC ransac(ransac_model);
	ransac.setDistanceThreshold(25);
	ransac.computeModel();
	ransac.refineModel();

	Eigen::VectorXf coeffs;
	ransac.getModelCoefficients(coeffs);

	double scale = ransac_model->m_scale;
	Eigen::Matrix3f rot = ransac_model->m_rot;
	Eigen::Vector3f tran = ransac_model->m_tran;

	Eigen::Matrix4f Tresult;
	Tresult.row(0).matrix() = coeffs.segment<4>(0);
	Tresult.row(1).matrix() = coeffs.segment<4>(4);
	Tresult.row(2).matrix() = coeffs.segment<4>(8);
	Tresult.row(3).matrix() = coeffs.segment<4>(12);
	
	//TransformPointCloud(vCam_src_, vCam_tgt_, Tresult);

	//modify EO
	double dePhi = std::atan2(rot(2,1), rot(2,2));
	double deOmega = std::asin(rot(2,0));
	double deKappa = std::atan2(rot(1,0), rot(0,0));


	ofstream oEos;
	char* oeo = argv[4];
	oEos.open("EO.txt");
	oEos << "phi omega kappa tx ty tz f" << std::endl;
	for (size_t i = 0; i < cam_num; i++)
	{
		oEos << fixed << setprecision(6) << *(eos + 7 * i) + dePhi << " " << *(eos + 7 * i + 1) + deOmega << " " << *(eos + 7 * i + 2) + deKappa << " "
			<< vCam_tgt_[i].x << " " << vCam_tgt_[i].y << " " << vCam_tgt_[i].z << " " << *(eos + 7 * i + 6) << std::endl;
		//*(eos + 7 * i) += dePhi;
		//*(eos + 7 * i + 1) += deOmega;
		//*(eos + 7 * i + 2) += deKappa;
		//*(eos + 7 * i + 3) += tran(0);
		//*(eos + 7 * i + 4) += tran(1);
		//*(eos + 7 * i + 5) += tran(2);
	}
	oEos.close();


	float *pts;
	int *rgb;
	float *norm;
	float *z_ptr;
	if (!is_dense)
	{
		rgb = new int[pts_num * 3];
		pts = new float[pts_num * 2];
		z_ptr = new float[pts_num];
		for (size_t i = 0; i < pts_num; i++)
		{
			float r, g, b;
			//read pointcloud for dem generation
			iBundle >> *(pts + i * 2) >> *(pts + i * 2 + 1) >> *(z_ptr + i) >> r >> g >> b;
			std::getline(iBundle, line);
			std::getline(iBundle, line);
			ptXYZ pt;

			//read pointcloud for coordinate transformation from model space to UTM coordniate system
			pt.x = *(pts + i * 2);
			pt.y = *(pts + i * 2 + 1);
			pt.z = *(z_ptr + i);

			*(rgb + i * 3) = static_cast<int>(r);
			*(rgb + i * 3 + 1) = static_cast<int>(g);
			*(rgb + i * 3 + 2) = static_cast<int>(b);
			pts_src.push_back(pt);
		}
	}
	else
	{
		char * densecloud = argv[6];
		ifstream iDense;
		iDense.open(densecloud);

		std::getline(iDense, line);
		if (line.find_first_of('p') != 0)
		{
			std::cout << "wrong ply file: " << densecloud << std::endl;
			return -1;
		}
		getline(iDense, line);

		getline(iDense, line);
		char *p = (char*)line.data();
		strtok(p, " ");
		strtok(NULL, " ");
		pts_num = atoi((const char *)strtok(NULL, " "));

		for (size_t i = 0; i < 10; i++)
		{
			std::getline(iDense, line);
		}
		
		if (line.find_first_of('e') != 0)
		{
			std::cout << "wrong ply head format: " << densecloud << std::endl;
			return -1;
		}
		
		pts = new float[pts_num * 2];
		rgb = new int[pts_num * 3];
		norm = new float[pts_num * 3];
		z_ptr = new float[pts_num];
		for (size_t i = 0; i < pts_num; i++)
		{
			float nx, ny, nz;
			float r, g, b;
			//read pointcloud for dem generation
			iDense >> *(pts + i * 2) >> *(pts + i * 2 + 1) >> *(z_ptr + i)
				>> *(norm + i * 3) >> *(norm + i * 3 + 1) >> *(norm + i * 3 + 2)
				>> r >> g >> b;

			*(rgb + i * 3) = static_cast<int>(r);
			*(rgb + i * 3 + 1) = static_cast<int>(g);
			*(rgb + i * 3 + 2) = static_cast<int>(b);
			ptXYZ pt;

			//read pointcloud for coordinate transformation from model space to UTM coordniate system
			pt.x = *(pts + i * 2);
			pt.y = *(pts + i * 2 + 1);
			pt.z = *(z_ptr + i);
			pts_src.push_back(pt);
		}

		iDense.close();
	}
	iBundle.close();


	//transform coordinate from model to UTM
	TransformPointCloud(pts_src, pts_tgt, Tresult);
	for (size_t i = 0; i < pts_tgt.size(); i++)
	{
		if (minx > pts_tgt[i].x) 
			minx = pts_tgt[i].x;
		if (miny > pts_tgt[i].y) 
			miny = pts_tgt[i].y;
		if (maxx < pts_tgt[i].x) 
			maxx = pts_tgt[i].x;
		if (maxy < pts_tgt[i].y) 
			maxy = pts_tgt[i].y;
	}

	//output the pointcloud in UTM
	std::string pointcloud_tgt = argv[5];
	ofstream ofile;
	ofile.open(pointcloud_tgt);
	if (!ofile)
	{
		std::cout << "Cannot write pointcloud file: " << pointcloud_tgt << std::endl;
		return -1;
	}
	ofile << "ply\nformat ascii 1.0\nelement vertex " << pts_num <<"\nproperty float x\nproperty float y\nproperty float z\n"
		<< "property float nx\nproperty float ny\nproperty float nz\n" << " property uchar diffuse_red\nproperty uchar diffuse_green\nproperty uchar diffuse_blue\nend_header\n";
	float *pts_tgt_ = new float[2 * pts_num];
	float *z_tgt = new float[pts_num];
	
	for (size_t i = 0; i < pts_num; i++)
	{
		*(pts_tgt_ + 2 * i) = pts_tgt[i].x;
		*(pts_tgt_ + 2 * i + 1) = pts_tgt[i].y;
		*(z_tgt + i) = pts_tgt[i].z;
		ofile << fixed << setprecision(6) << pts_tgt[i].x + deltaX << " " << pts_tgt[i].y +deltaY << " " << pts_tgt[i].z ;
		if(!is_dense)
			ofile << " 0.0 0.0 0.0 ";
		else
			ofile << fixed << setprecision(3) << " " << *(norm+i*3) << " " << *(norm + i *3 +1) << " " << *(norm + i *3 +2) << " ";
		ofile << (int)*(rgb+i*3) << " " << (int)*(rgb + i *3 +1) << " " << (int)*(rgb + i *3 +2) << std::endl;
	}
	ofile.close();

	//generate dem
	flann::Matrix<float> features(pts_tgt_, pts_num, 2);
	flann::Index<flann::L2<float> > index(features, flann::KDTreeSingleIndexParams());
	index.buildIndex();
	int minX = static_cast<int>(minx);
	int minY = static_cast<int>(miny);
	int width = static_cast<int>((maxx -  minx) / 0.1);
	int height = static_cast<int>((maxy - miny) / 0.1);

	char * dem_tgt = argv[3];
	ofile.open(dem_tgt);
	if (!ofile)
	{
		std::cout << "Cannot write dem file: " << dem_tgt << std::endl;
		return -1;
	}

	ofile << fixed << setprecision(6) << minx+deltaX << " " << miny+deltaY << " 0.1 0.1 " << width << " " << height 
		<<" " << deltaX << " " << deltaY << " " << zone <<std::endl;
	double radius = 2;
	for (size_t i = 0; i < width; i++)
	{
		float x = minX + i * 0.1;

		for (size_t j = 0; j < height; j++)
		{
			float y = minY + j * 0.1;
			std::vector< std::vector<int> > indices;
			std::vector< std::vector<float> > dists;
			float *pos = new float[2];
			*pos = x;
			*(pos+1) = y;
			flann::Matrix<float> query(pos, 1, 2);
			int count = index.radiusSearch(query, indices, dists, radius, flann::SearchParams(10));
			double t_pz = 0;
			double t_p = 0;
			if (count > 6)
			{
				double z;
				bool z_flag = false;
				for (size_t k = 0; k < count; k++)
				{
					double dis = dists[0][k];
					if (dis == 0)
					{
						z = *(z_tgt + indices[0][0]);
						ofile << fixed << setprecision(6) << z << " ";
						z_flag = true;
						break;
					}
					double pp = 1. / dis;
					t_pz += *(z_tgt + indices[0][k]) / dis;
					t_p += pp;
				}
				if (!z_flag)
				{
					z = t_pz / t_p;
					ofile << z << " ";
				}
			}
			else
			{
				ofile << "-9999.999999 ";
			}
		}
		ofile << std::endl;
	}
	ofile.close();
	return 0;
}
