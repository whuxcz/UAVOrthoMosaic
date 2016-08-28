#include "stdio.h"
#include "math.h"
#include "malloc.h"
#include "matrix.h"
#include <fstream>
#include <sstream>
#include <stdlib.h> 
#include <iostream>       
#include <string>         
#include <cstddef>         
#include <vector>
#include "imgutil.hpp"


#define grid_size 0.1
#define PI 3.1415926
#define INF -9999.999999

#define _MAX_PATH 260
#define _MAX_DRIVE  3
#define _MAX_DIR 256
#define _MAX_FNAME 256
#define _MAX_EXT 256


struct DEMPoint{
	double X;
	double Y;
	double Z;
};

struct PixelPoint
{
	double x;
	double y;
};

DEMPoint *dem_pt;
PixelPoint *px_pt;
int rows, cols;

int gridminx, gridmaxx, gridminy, gridmaxy;
double UTM_X0, UTM_Y0, offside_X, offside_Y;
int X0, Y0, ZONE;
double dX, dY;
double aver;
double xy[2];

double l_X0, l_Y0;
int i_count, j_count;

void GetDEM(const char *DEM_PathName)
{
	std::ifstream fin(DEM_PathName, std::ios::in);
	fin >> UTM_X0 >> UTM_Y0 >> dX >> dY >> rows >> cols >> offside_X >> offside_Y >> ZONE;
	X0 = int(UTM_X0 - offside_X);
	Y0 = int(UTM_Y0 - offside_Y);

	int num = rows*cols;
	dem_pt = new DEMPoint[num];
	double sum_z = 0;
	int z_valid = 0;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{

			dem_pt[i*cols + j].X = X0 + i*dX;
			dem_pt[i*cols + j].Y = Y0 + j*dY;
			if (!fin.eof())
			{
				fin >> dem_pt[i*cols + j].Z;
				if (dem_pt[i*cols + j].Z != INF)
				{
					sum_z += dem_pt[i*cols + j].Z;
					z_valid++;
				}
			}


		}
	}
	aver = sum_z / z_valid;

	fin.close();
}
void Interpolation(const char *INPUT_PathName, const char *DOM_PathName)
{
	IplImage* inImg;
	IplImage* domImg;

	inImg = iu_LoadJpgImage(INPUT_PathName);



	domImg = iu_CreateImage(j_count, i_count, inImg->nChannels);
	IuScalar temp;
	temp.val[0] = 0; temp.val[1] = 0; temp.val[2] = 0;
	for (int i = 0; i < domImg->height; i++)
	{
		for (int j = 0; j < domImg->width; j++)
		{
			iu_Set2D(domImg, i, j, temp);

		}

	}




	int int_x, int_y;
	double wx, wy;
	IuScalar s1, s2, s3, s4;
	IuScalar result;


	double  R, G, B;
	for (int i = 0; i < inImg->height; i++)   
	{
		for (int j = 0; j < inImg->width; j++) 
		{
			
			R = 0;
			G = 0;
			B = 0;

			int_x = int(px_pt[i*inImg->width + j].x);

			int_y = int(px_pt[i*inImg->width + j].y);
			if (int_x >= 0 && int_x < inImg->width - 1 && int_y >= 0 && int_y < inImg->height - 1)

			{
				s1 = iu_Get2D(inImg, inImg->height - 1 - int_y, int_x);
				s2 = iu_Get2D(inImg, inImg->height - 1 - int_y, int_x + 1);
				s3 = iu_Get2D(inImg, inImg->height - 1 - int_y - 1, int_x);
				s4 = iu_Get2D(inImg, inImg->height - 1 - int_y - 1, int_x + 1);

				wx = px_pt[i*inImg->width + j].x - int_x;
				wy = px_pt[i*inImg->width + j].y - int_y;


				B = (1 - wx) * (s1.val[0] * (1 - wy) + s2.val[0] * wy) + wx * (s3.val[0] * (1 - wy) + s4.val[0] * wy);
				G = (1 - wx) * (s1.val[1] * (1 - wy) + s2.val[1] * wy) + wx * (s3.val[1] * (1 - wy) + s4.val[1] * wy);
				R = (1 - wx) * (s1.val[2] * (1 - wy) + s2.val[2] * wy) + wx * (s3.val[2] * (1 - wy) + s4.val[2] * wy);



				result.val[0] = B;
				result.val[1] = G;
				result.val[2] = R;
				if (i < i_count - 1 && j < j_count - 1)
					iu_Set2D(domImg, i, j, result);

			}


		}

	}

	double tiepoints[6] = { 0.0, 0.0, 0.0, l_X0 + offside_X, l_Y0 + offside_Y, 0.0 };
	double pixscale[3] = { xy[0], xy[1], 0 };
	
	IplImage* pRgbaImg = iu_CreateImage(domImg->width, domImg->height, 4);
	for (int i = 0; i < domImg->height; i++)
	{
		for (int j = 0; j < domImg->width; j++)
		{
			IuScalar s = iu_Get2D(domImg, i, j);
			if (s.val[0] == 0 && s.val[1] == 0 && s.val[2] == 0)
			{
				s.val[3] = 0;
			}
			else
			{
				s.val[3] = 255;
			}
			iu_Set2D(pRgbaImg, i, j, s);
		}
	}

	iu_SaveGTifImage(DOM_PathName, pRgbaImg, tiepoints, pixscale, ZONE);
	delete[](inImg->imageData);
	delete inImg;
	delete[](domImg->imageData);
	delete domImg;
	delete[](pRgbaImg->imageData);
	delete pRgbaImg;


	


}
double Deminterpolation(double x, double y, const char *DEM_PathName)
{

	DEMPoint  trdem_pt[9], rdem_pt[4];
	int i, j, k, surround = 0;

	int w_rows = int((x - X0) / grid_size);
	int  w_cols = int((y - Y0) / grid_size);



	if ((w_cols >= 0 && w_cols < (cols - 1)) && (w_rows >= 0 && w_rows < (rows - 1)))
	{

		for (i = 0; i < 9; i++)
		{
			trdem_pt[i].Z = INF;

		}

		for (i = 0; i < 4; i++)
		{
			rdem_pt[i].Z = INF;
		}

		for (i = -1; i < 2 && (w_rows + i) >= 0 && (w_rows + i) < rows; i++)
		{
			for (j = -1; j < 2 && (w_cols + j) >= 0 && (w_cols + j) < cols; j++)
			{
				trdem_pt[(i + 1) * 3 + j + 1].X = dem_pt[(w_rows + i)* cols + w_cols + j].X;
				trdem_pt[(i + 1) * 3 + j + 1].Y = dem_pt[(w_rows + i)* cols + w_cols + j].Y;
				trdem_pt[(i + 1) * 3 + j + 1].Z = dem_pt[(w_rows + i)* cols + w_cols + j].Z;
				surround++;


			}
		}
		k = 0;
		if (trdem_pt[4].Z == INF&&trdem_pt[5].Z == INF&&trdem_pt[7].Z == INF&&trdem_pt[8].Z == INF)return 0;

		else
		{
			if (trdem_pt[4].Z != INF)
			{
				rdem_pt[k].X = trdem_pt[4].X;
				rdem_pt[k].Y = trdem_pt[4].Y;
				rdem_pt[k].Z = trdem_pt[4].Z;
				k++;
			}
			if (trdem_pt[5].Z != INF)
			{
				rdem_pt[k].X = trdem_pt[5].X;
				rdem_pt[k].Y = trdem_pt[5].Y;
				rdem_pt[k].Z = trdem_pt[5].Z;
				k++;
			}
			if (trdem_pt[7].Z != INF)
			{
				rdem_pt[k].X = trdem_pt[7].X;
				rdem_pt[k].Y = trdem_pt[7].Y;
				rdem_pt[k].Z = trdem_pt[7].Z;
				k++;
			}
			if (trdem_pt[8].Z != INF)
			{
				rdem_pt[k].X = trdem_pt[8].X;
				rdem_pt[k].Y = trdem_pt[8].Y;
				rdem_pt[k].Z = trdem_pt[8].Z;
				k++;
			}



		}
		if (surround >= 4)
		{
			for (i = 1, j = k; i < 9 && i != 4 && i != 5 && i != 7 && i != 8 && j < 4 && trdem_pt[i].Z != INF; i++)
			{
				rdem_pt[j].X = trdem_pt[i].X;
				rdem_pt[j].Y = trdem_pt[i].Y;
				rdem_pt[j].Z = trdem_pt[i].Z;
				j++;
			}

		}


		if (rdem_pt[0].Z == INF || rdem_pt[1].Z == INF || rdem_pt[2].Z == INF || rdem_pt[3].Z == INF)return 0;

		else
		{
			double m[4][4], mt[4][4], mtm[4][4], mtmn[4][4], mtz[4][1], xx[4][1], z[4][1];

			for (i = 0; i < 4; i++)
			{
				m[i][0] = 1;
				m[i][1] = rdem_pt[i].X;
				m[i][2] = rdem_pt[i].Y;
				m[i][3] = rdem_pt[i].X*rdem_pt[i].Y;
				z[i][0] = rdem_pt[i].Z;
			}
			for (i = 0; i < 4; i++)
			{
				for (j = 0; j < 4; j++)
				{
					mt[i][j] = m[j][i];
				}
			}
			for (i = 0; i < 4; i++)
			{
				for (j = 0; j < 4; j++)
				{
					mtm[i][j] = 0;
					for (k = 0; k < 4; k++)
					{
						mtm[i][j] += mt[i][k] * m[k][j];
					}
				}
			}

			Matrix mtmn_1(4, 4), mtm_1(4, 4);

			for (i = 0; i < 4; i++)
			{
				for (j = 0; j < 4; j++)
				{
					mtm_1.change(i, j, mtm[i][j]);
				}

			}
			mtmn_1 = mtm_1.inverse();
			for (i = 0; i < 4; i++)
			{
				for (j = 0; j < 4; j++)
				{
					mtmn[i][j] = mtmn_1.get(i, j);
				}

			}
			for (i = 0; i < 4; i++)
			{
				for (j = 0; j < 1; j++)
				{
					mtz[i][j] = 0;
					for (k = 0; k < 4; k++)
					{
						mtz[i][j] += mt[i][k] * z[k][j];
					}
				}
			}
			for (i = 0; i < 4; i++)
			{
				for (j = 0; j < 1; j++)
				{
					xx[i][j] = 0;
					for (k = 0; k < 4; k++)
					{
						xx[i][j] += mtmn[i][k] * mtz[k][j];
					}
				}
			}


			return  xx[0][0] + xx[1][0] * x + xx[2][0] * y + xx[3][0] * x*y;

		}

	}
	else
	{
		return 0;
	}

	free(trdem_pt);
	free(rdem_pt);



}
void Sample_Cal(const char *path[], const char *m_fileop, double xy[], const char *DEM_PathName, int size)
{
	GetDEM(DEM_PathName);
	double Xs, Ys, Zs, phi, omega, kappa, x0, y0, f, l_X0, l_Y0, l_X1, l_Y1;

	
	std::ifstream fin(m_fileop, std::ios::in);
	std::string info;
	getline(fin, info);


	double(*ao)[7];
	ao = (double(*)[7])malloc(size * 7 * sizeof(double));
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < 7; j++)
		{

			fin >> ao[i][j];
		}
	}
	xy[0] = 0; xy[1] = 0;
	for (size_t t = 0; t < size; t++)
	{
		IplImage* inImg;
		inImg = iu_LoadJpgImage(path[t]);
		phi = ao[t][0], omega = ao[t][1], kappa = ao[t][2], Xs = ao[t][3], Ys = ao[t][4], Zs = ao[t][5], x0 = 2000, y0 = 1500, f = ao[t][6];
		double a1, a2, a3, b1, b2, b3, c1, c2, c3;
		a1 = cos(phi)*cos(kappa) - sin(phi)*sin(omega)*sin(kappa);
		a2 = -cos(phi)*sin(kappa) - sin(phi)*sin(omega)*cos(kappa);
		a3 = -sin(phi)*cos(omega);
		b1 = cos(omega)*sin(kappa);
		b2 = cos(omega)*cos(kappa);
		b3 = -sin(omega);
		c1 = sin(phi)*cos(kappa) + cos(phi)*sin(omega)*sin(kappa);
		c2 = -sin(phi)*sin(kappa) + cos(phi)*sin(omega)*cos(kappa);
		c3 = cos(phi)*cos(omega);

		double L, L1, L2;
		L = -(a3*Xs + b3*Ys + c3*Zs);
		L1 = -(a1*Xs + b1*Ys + c1*Zs);
		L2 = -(a2*Xs + b2*Ys + c2*Zs);
		double Z0 = aver, Z1 = aver + 0.2;
		for (int k = 0; (k < 100) && (fabs(Z1 - Z0)>0.01); k++)
		{
			Z0 = Z1;

			l_X0 = (Z0 - Zs)*(-inImg->width / 2 * a1 + (inImg->height / 2) * a2 - a3*f) / (-inImg->width / 2 * c1 + (inImg->height / 2) * c2 - f*c3) + Xs;
			l_Y0 = (Z0 - Zs)*(-inImg->width / 2 * b1 + (inImg->height / 2) * b2 - b3*f) / (-inImg->width / 2 * c1 + (inImg->height / 2) * c2 - f*c3) + Ys;

			Z1 = Deminterpolation(l_X0, l_Y0, DEM_PathName);
			if (Z1 == 0)break;


		}
		Z0 = aver, Z1 = aver + 0.2;
		for (int k = 0; (k < 100) && (fabs(Z1 - Z0)>0.01); k++)
		{
			Z0 = Z1;

			l_X1 = (Z0 - Zs)*(inImg->width / 2 * a1 + (-inImg->height / 2) * a2 - a3*f) / (inImg->width / 2 * c1 + (-inImg->height / 2) * c2 - f*c3) + Xs;
			l_Y1 = (Z0 - Zs)*(inImg->width / 2 * b1 + (-inImg->height / 2) * b2 - b3*f) / (inImg->width / 2 * c1 + (-inImg->height / 2) * c2 - f*c3) + Ys;

			Z1 = Deminterpolation(l_X1, l_Y1, DEM_PathName);
			if (Z1 == 0)break;


		}





		double m_xTemp = fabs(cos(kappa)*(l_X1 - l_X0) + sin(kappa)*(l_Y1 - l_Y0)) / inImg->width;
		double m_yTemp = fabs(-sin(kappa)*(l_X1 - l_X0) + cos(kappa)*(l_Y1 - l_Y0)) / inImg->height;
		delete[](inImg->imageData);
		delete inImg;
		if (m_xTemp > xy[0])xy[0] = m_xTemp;
		if (m_yTemp > xy[1])xy[1] = m_yTemp;
	}
	double aa = xy[0];
	double bb = xy[1];
	free(ao);

}

void DomcorrectionSets(const char *path[], const char *pathout[], const char *m_fileop, const char *DEM_PathName, const char *DOM_info, int size)
{


	Sample_Cal(path, m_fileop, xy, DEM_PathName, size);

	std::ifstream fin(m_fileop, std::ios::in);

	


	IplImage* inImg;

	std::string info1;
	std::getline(fin, info1);


	double(*ao)[7];
	ao = (double(*)[7])malloc(size * 7 * sizeof(double));

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < 7; j++)
		{

			fin >> ao[i][j];
		}
	}


	std::ofstream out(DOM_info, std::ios::app);
	out << "X0  Y0  X_Resolution Y_Resolution  A" << endl;
	for (int t = 0; t < size; t++)
	{

		inImg = iu_LoadJpgImage(path[t]);
		double Xs, Ys, Zs, phi, omega, kappa, x0, y0, f, l_X1, l_Y1;
		phi = ao[t][0], omega = ao[t][1], kappa = ao[t][2], Xs = ao[t][3], Ys = ao[t][4], Zs = ao[t][5], x0 = 2000, y0 = 1500, f = ao[t][6];
		int r_rows = inImg->height;
		int r_cols = inImg->width;
		px_pt = new PixelPoint[r_rows*r_cols];




		double a1, a2, a3, b1, b2, b3, c1, c2, c3;
		a1 = cos(phi)*cos(kappa) - sin(phi)*sin(omega)*sin(kappa);
		a2 = -cos(phi)*sin(kappa) - sin(phi)*sin(omega)*cos(kappa);
		a3 = -sin(phi)*cos(omega);
		b1 = cos(omega)*sin(kappa);
		b2 = cos(omega)*cos(kappa);
		b3 = -sin(omega);
		c1 = sin(phi)*cos(kappa) + cos(phi)*sin(omega)*sin(kappa);
		c2 = -sin(phi)*sin(kappa) + cos(phi)*sin(omega)*cos(kappa);
		c3 = cos(phi)*cos(omega);

		double L, L1, L2;
		L = -(a3*Xs + b3*Ys + c3*Zs);
		L1 = -(a1*Xs + b1*Ys + c1*Zs);
		L2 = -(a2*Xs + b2*Ys + c2*Zs);
		double Z0 = aver, Z1 = aver + 0.2;
		for (int k = 0; (k < 100) && (fabs(Z1 - Z0)>0.01); k++)
		{
			Z0 = Z1;

			l_X0 = (Z0 - Zs)*(-inImg->width / 2 * a1 + (inImg->height / 2) * a2 - a3*f) / (-inImg->width / 2 * c1 + (inImg->height / 2) * c2 - f*c3) + Xs;
			l_Y0 = (Z0 - Zs)*(-inImg->width / 2 * b1 + (inImg->height / 2) * b2 - b3*f) / (-inImg->width / 2 * c1 + (inImg->height / 2) * c2 - f*c3) + Ys;

			Z1 = Deminterpolation(l_X0, l_Y0, DEM_PathName);
			if (Z1 == 0)break;


		}
		Z0 = aver, Z1 = aver + 0.2;
		for (int k = 0; (k < 100) && (fabs(Z1 - Z0)>0.01); k++)
		{
			Z0 = Z1;

			l_X1 = (Z0 - Zs)*(inImg->width / 2 * a1 + (-inImg->height / 2) * a2 - a3*f) / (inImg->width / 2 * c1 + (-inImg->height / 2) * c2 - f*c3) + Xs;
			l_Y1 = (Z0 - Zs)*(inImg->width / 2 * b1 + (-inImg->height / 2) * b2 - b3*f) / (inImg->width / 2 * c1 + (-inImg->height / 2) * c2 - f*c3) + Ys;

			Z1 = Deminterpolation(l_X1, l_Y1, DEM_PathName);
			if (Z1 == 0)break;


		}

		//cvReleaseImage(&inImg);
		delete[](inImg->imageData);
		delete inImg;
		double l_x, l_y, l_z;


		j_count = int(fabs(cos(kappa)*(l_X1 - l_X0) + sin(kappa)*(l_Y1 - l_Y0)) / xy[0]);
		i_count = int(fabs(-sin(kappa)*(l_X1 - l_X0) + cos(kappa)*(l_Y1 - l_Y0)) / xy[1]);


		out << l_X0 << " " << l_Y0 << " " << xy[0] << " " << xy[1] << " " << kappa << endl;


		for (int i = 0; i < i_count; i++)
		{
			for (int j = 0; j < j_count; j++)
			{

				l_x = l_X0 + xy[0] * j*cos(kappa) - (-xy[1] * i)*sin(kappa);
				l_y = l_Y0 + xy[0] * j*sin(kappa) + (-xy[1] * i)*cos(kappa);

				if ((l_x > X0 && l_x<(X0 + dX*rows)) && (l_y>Y0 && l_y < (Y0 + dY*cols)))
				{
					l_z = Deminterpolation(l_x, l_y, DEM_PathName);

					if (l_z == 0)
					{
						px_pt[i* r_cols + j].x = -1;
						px_pt[i* r_cols + j].y = -1;
					}
					else
					{
						px_pt[i* r_cols + j].x = -f*(l_x*a1 + l_y*b1 + l_z*c1 + L1) / (l_x*a3 + l_y*b3 + l_z*c3 + L) + x0;
						px_pt[i* r_cols + j].y = -f*(l_x*a2 + l_y*b2 + l_z*c2 + L2) / (l_x*a3 + l_y*b3 + l_z*c3 + L) + y0;
					}
				}
				else
				{
					px_pt[i* r_cols + j].x = -1;
					px_pt[i* r_cols + j].y = -1;
				}


			}

		}



		Interpolation(path[t], pathout[t]);
		free(px_pt);

	}

}



void _split_whole_name(const char *whole_name, char *fname, char *ext)
{
	char *p_ext;
	char *p = new char[strlen(whole_name) + 1];
	strcpy(p, whole_name);
	p_ext = rindex(p, '.');
	if (NULL != p_ext)
	{
		strcpy(ext, p_ext);
		snprintf(fname, p_ext - whole_name + 1, "%s", whole_name);
	}
	else
	{
		ext[0] = '\0';
		strcpy(fname, whole_name);
	}
}

void _splitpath(const char *path, char *drive, char *dir, char *fname, char *ext)
{
	char *p_whole_name;

	drive[0] = '\0';
	if (NULL == path)
	{
		dir[0] = '\0';
		fname[0] = '\0';
		ext[0] = '\0';
		return;
	}

	if ('/' == path[strlen(path)])
	{
		strcpy(dir, path);
		fname[0] = '\0';
		ext[0] = '\0';
		return;
	}
	char *p = new char[strlen(path) + 1];
	strcpy(p, path);
	p_whole_name = rindex(p, '/');
	if (NULL != p_whole_name)
	{
		p_whole_name++;
		_split_whole_name(p_whole_name, fname, ext);

		snprintf(dir, p_whole_name - path, "%s", path);
	}
	else
	{
		_split_whole_name(path, fname, ext);
		dir[0] = '\0';
	}
}




int main(int argc, char* argv[])
{	
	if(argc != 7)
	{
		std::cout << "[Usage] MutiDOMgenerate: InputImgList outputImgSets InputEO InputDem_file OutputDominfo InputImgList_size " << std::endl;
		return -1;
	}
	const char* input;
	const char* output;
	const char* ao_input;
	const char* dem_input;
	const char* dominfo;
	const char* size;
	

	input = argv[1];
	output = argv[2];
	ao_input = argv[3];
	dem_input = argv[4];
	dominfo = argv[5];
	size = argv[6];
	

	int _size;
	
	_size = atoi(size);

	
	
	std::string* path=new std::string[_size];
	const char** path_input = new const char*[_size];
	const char** path_out1 = new const char*[_size];
	std::vector<std::string> path_out;
	std::ifstream fin(input, std::ios::in);
	for (size_t i = 0; i < _size ; i++)
	{
		getline(fin,path[i]);
		
	}




	char full_path[_MAX_PATH];
	char drive[_MAX_DRIVE];
	char dir[_MAX_DIR];
	char fname[_MAX_FNAME];
	char ext[_MAX_EXT];

	std::string dom = "-dom";
	std::string temp_dom;
	const char*  _fname;
  
	for (size_t i = 0; i < _size; i++)
	{
		path_input[i] = path[i].c_str();
		std::string path_output = path[i];
		path_output.erase(path_output.begin()+ path_output.find_last_of("."), path_output.end());
		path_output += "-dom.tif";
		//_splitpath(path_input[i], drive, dir, fname, ext);

		//std::stringextt = ".tif";
		//std::stringend = string(fname) + dom + extt;

		//std::stringheader = output;

		//std::stringpath_output = header + end;
		path_out.push_back(path_output);

		path_out1[i] = path_out[i].data();
		
	}
	DomcorrectionSets(path_input, path_out1, ao_input, dem_input, dominfo, _size);

	return 0;
}
