#ifndef __keys2a_flann_h__
#define __keys2a_flann_h__

#include <vector>
#include <stdio.h>

#include <zlib.h>
#include <algorithm>

#include <flann/flann.hpp>

using namespace flann;

typedef struct 
{	
	float x;
	float y;
	float h;
	float r;
}photoPosition;

class Keypoint {
public:
	Keypoint(float x, float y, float scale, float ori, short int *d) :
		m_x(x), m_y(y), m_scale(scale), m_ori(ori), m_d(d)
	{ }

	float m_x, m_y;        /* Subpixel location of keypoint. */
	float m_scale, m_ori;  /* Scale and orientation (range [-PI,PI]) */
	short int *m_d;    /* Vector of descriptor values */
};

/* Data struct for matches */
class KeypointMatch {
public:
	KeypointMatch()
	{ }

#if 0
	KeypointMatch(int idx1, int idx2, float x1, float y1, float x2, float y2) :
		m_idx1(idx1), m_idx2(idx2), m_x1(x1), m_y1(y1), m_x2(x2), m_y2(y2)
	{ }
#endif

	KeypointMatch(int idx1, int idx2) :
		m_idx1(idx1), m_idx2(idx2)
	{ }

	int m_idx1, m_idx2;
	// float m_x1, m_y1;
	// float m_x2, m_y2;
};


typedef struct {
	float x, y;
	float scale;
	float orient;
} keypt_t;

/* Returns the number of keys in a file */
int GetNumberOfKeys(const char *filename);

/* This reads a keypoint file from a given filename and returns the list
* of keypoints. */
int ReadKeyFile(const char *filename, unsigned char **keys,	keypt_t **info = NULL);

int ReadKeyPositions(const char *filename, keypt_t **info);

/* Read keypoints from the given file pointer and return the list of
* keypoints.  The file format starts with 2 integers giving the total
* number of keypoints and the size of descriptor vector for each
* keypoint (currently assumed to be 128). Then each keypoint is
* specified by 4 floating point numbers giving subpixel row and
* column location, scale, and orientation (in radians from -PI to
* PI).  Then the descriptor vector for each keypoint is given as a
* list of integers in range [0,255]. */
int ReadKeys(FILE *fp, unsigned char **keys, keypt_t **info = NULL);
int ReadKeysGzip(gzFile fp, unsigned char **keys, keypt_t **info = NULL);


/* Read keys using MMAP to speed things up */
std::vector<Keypoint *> ReadKeysMMAP(FILE *fp);

int ReadPosition(const char *filename, std::vector<photoPosition> &pos, float f_mm, float ccd_width);


/* Create a search tree for the given set of keypoints based on flann*/
Index<L2<unsigned char> > *CreateSearchTree(int num_keys, unsigned char *keys);


/* Compute likely matches between two sets of keypoints */
std::vector<KeypointMatch> MatchKeys(int num_key1, unsigned char *k1,
	Index<L2<unsigned char> > index,
	double ratio = 0.6, int max_pts_visit=15);

#endif /*__keys2a_flann_h__*/
