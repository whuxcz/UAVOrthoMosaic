#include "keys2a_flann.h"
#include <vector>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <zlib.h>


int GetNumberOfKeysNormal(FILE *fp)
{
	int num, len;

	if (fscanf(fp, "%d %d", &num, &len) != 2) {
		printf("Invalid keypoint file.\n");
		return 0;
	}

	return num;
}

int GetNumberOfKeysGzip(gzFile fp)
{
	int num, len;

	char header[256];
	gzgets(fp, header, 256);

	if (sscanf(header, "%d %d", &num, &len) != 2) {
		printf("Invalid keypoint file.\n");
		return 0;
	}

	return num;
}

/* Returns the number of keys in a file */
int GetNumberOfKeys(const char *filename)
{
	FILE *file;

	file = fopen(filename, "r");
	if (!file) {
		/* Try to file a gzipped keyfile */
		char buf[1024];
		sprintf(buf, "%s.gz", filename);
		gzFile gzf = gzopen(buf, "rb");

		if (gzf == NULL) {
			printf("Could not open file: %s\n", filename);
			return 0;
		}
		else {
			int n = GetNumberOfKeysGzip(gzf);
			gzclose(gzf);
			return n;
		}
	}

	int n = GetNumberOfKeysNormal(file);
	fclose(file);
	return n;
}

int ReadPosition(const char *filename, std::vector<photoPosition> &pos, float f_mm, float ccd_width)
{	
	FILE *file;

	int num;
	float x,y,h,r;

	file = fopen(filename, "r");
	if(!file)
	{	
		printf("Cannot open postion file: %s\n", filename);
	}	

	if (fscanf(file, "%d", &num) != 1) {
		printf("Invalid position file\n");
		return 0;
	}
	
	for(size_t i = 0; i < num; i++)
	{
		if (fscanf(file, "%f %f %f", &x, &y, &h) != 3) 
		{
			printf("Invalid position file\n");
			return 0;
		}
		r = ccd_width * h / f_mm;
		photoPosition position;
		position.x = x;
		position.y = y;
		position.h = h;
		position.r = r;
		pos.push_back(position);
	}
	return num;
}

/* This reads a keypoint file from a given filename and returns the list
* of keypoints. */
int ReadKeyFile(const char *filename, unsigned char **keys, keypt_t **info)
{
	FILE *file;

	file = fopen(filename, "r");
	if (!file) {
		/* Try to file a gzipped keyfile */
		char buf[1024];
		sprintf(buf, "%s.gz", filename);
		gzFile gzf = gzopen(buf, "rb");

		if (gzf == NULL) {
			printf("Could not open file: %s\n", filename);
			return 0;
		}
		else {
			int n = ReadKeysGzip(gzf, keys, info);
			gzclose(gzf);
			return n;
		}
	}

	int n = ReadKeys(file, keys, info);
	fclose(file);
	return n;
	// return ReadKeysMMAP(file);
}

#if 0
/* Read keys using MMAP to speed things up */
std::vector<Keypoint *> ReadKeysMMAP(FILE *fp)
{
	int i, j, num, len, val, n;

	std::vector<Keypoint *> kps;

	struct stat sb;

	/* Stat the file */
	if (fstat(fileno(fp), &sb) < 0) {
		printf("[ReadKeysMMAP] Error: could not stat file\n");
		return kps;
	}

	char *file = (char *)mmap(NULL, sb.st_size, PROT_READ, MAP_SHARED,
		fileno(fp), 0);

	char *file_start = file;

	if (sscanf(file, "%d %d%n", &num, &len, &n) != 2) {
		printf("[ReadKeysMMAP] Invalid keypoint file beginning.");
		return kps;
	}

	file += n;

	if (len != 128) {
		printf("[ReadKeysMMAP] Keypoint descriptor length invalid "
			"(should be 128).");
		return kps;
	}

	for (i = 0; i < num; i++) {
		/* Allocate memory for the keypoint. */
		unsigned char *d = new unsigned char[len];
		float x, y, scale, ori;

		if (sscanf(file, "%f %f %f %f%n", &y, &x, &scale, &ori, &n) != 4) {
			printf("[ReadKeysMMAP] Invalid keypoint file format.");
			return kps;
		}

		file += n;

		for (j = 0; j < len; j++) {
			if (sscanf(file, "%d%n", &val, &n) != 1 || val < 0 || val > 255) {
				printf("[ReadKeysMMAP] Invalid keypoint file value.");
				return kps;
			}
			d[j] = (unsigned char)val;
			file += n;
		}

		kps.push_back(new Keypoint(x, y, scale, ori, d));
	}

	/* Unmap */
	if (munmap(file_start, sb.st_size) < 0) {
		printf("[ReadKeysMMAP] Error: could not unmap memory\n");
		return kps;
	}

	return kps;
}
#endif

/* Read keypoints from the given file pointer and return the list of
* keypoints.  The file format starts with 2 integers giving the total
* number of keypoints and the size of descriptor vector for each
* keypoint (currently assumed to be 128). Then each keypoint is
* specified by 4 floating point numbers giving subpixel row and
* column location, scale, and orientation (in radians from -PI to
* PI).  Then the descriptor vector for each keypoint is given as a
* list of integers in range [0,255]. */
int ReadKeys(FILE *fp, unsigned char **keys, keypt_t **info)
{
	int i, num, len;

	std::vector<Keypoint *> kps;

	if (fscanf(fp, "%d %d", &num, &len) != 2) {
		printf("Invalid keypoint file\n");
		return 0;
	}

	if (len != 128) {
		printf("Keypoint descriptor length invalid (should be 128).");
		return 0;
	}

	*keys = new unsigned char[128 * num + 8];

	if (info != NULL)
		*info = new keypt_t[num];

	unsigned char *p = *keys;
	for (i = 0; i < num; i++) {
		/* Allocate memory for the keypoint. */
		// short int *d = new short int[128];
		float x, y, scale, ori;

		if (fscanf(fp, "%f %f %f %f\n", &y, &x, &scale, &ori) != 4) {
			printf("Invalid keypoint file format.");
			return 0;
		}

		if (info != NULL) {
			(*info)[i].x = x;
			(*info)[i].y = y;
			(*info)[i].scale = scale;
			(*info)[i].orient = ori;
		}

		char buf[1024];
		for (int line = 0; line < 7; line++) {
			fgets(buf, 1024, fp);

			if (line < 6) {
				sscanf(buf,
					"%hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu "
					"%hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu",
					p + 0, p + 1, p + 2, p + 3, p + 4, p + 5, p + 6, p + 7, p + 8, p + 9,
					p + 10, p + 11, p + 12, p + 13, p + 14,
					p + 15, p + 16, p + 17, p + 18, p + 19);

				p += 20;
			}
			else {
				sscanf(buf,
					"%hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu",
					p + 0, p + 1, p + 2, p + 3, p + 4, p + 5, p + 6, p + 7);
				p += 8;
			}
		}
	}

	return num; // kps;
}

int ReadKeysGzip(gzFile fp, unsigned char **keys, keypt_t **info)
{
	int i, num, len;

	std::vector<Keypoint *> kps;
	char header[256];
	gzgets(fp, header, 256);

	if (sscanf(header, "%d %d", &num, &len) != 2) {
		printf("Invalid keypoint file.\n");
		return 0;
	}

	if (len != 128) {
		printf("Keypoint descriptor length invalid (should be 128).");
		return 0;
	}

	*keys = new unsigned char[128 * num + 8];

	if (info != NULL)
		*info = new keypt_t[num];

	unsigned char *p = *keys;
	for (i = 0; i < num; i++) {
		/* Allocate memory for the keypoint. */
		// short int *d = new short int[128];
		float x, y, scale, ori;
		char buf[1024];
		gzgets(fp, buf, 1024);

		if (sscanf(buf, "%f %f %f %f\n", &y, &x, &scale, &ori) != 4) {
			printf("Invalid keypoint file format.");
			return 0;
		}

		if (info != NULL) {
			(*info)[i].x = x;
			(*info)[i].y = y;
			(*info)[i].scale = scale;
			(*info)[i].orient = ori;
		}

		for (int line = 0; line < 7; line++) {
			char *str = gzgets(fp, buf, 1024);
			assert(str != Z_NULL);

			if (line < 6) {
				sscanf(buf,
					"%hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu "
					"%hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu",
					p + 0, p + 1, p + 2, p + 3, p + 4, p + 5, p + 6, p + 7, p + 8, p + 9,
					p + 10, p + 11, p + 12, p + 13, p + 14,
					p + 15, p + 16, p + 17, p + 18, p + 19);

				p += 20;
			}
			else {
				sscanf(buf,
					"%hhu %hhu %hhu %hhu %hhu %hhu %hhu %hhu",
					p + 0, p + 1, p + 2, p + 3, p + 4, p + 5, p + 6, p + 7);
				p += 8;
			}
		}
	}

	assert(p == *keys + 128 * num);

	return num; // kps;
}

Index<L2<unsigned char> > *CreateSearchTree(int num_keys, unsigned char *keys)
{
	Matrix<unsigned char> features(keys, num_keys, 128);
	Index<L2<unsigned char> > *index = new Index<L2<unsigned char> >(features, KDTreeIndexParams());
	index->buildIndex();
	return index;
}



std::vector<KeypointMatch> MatchKeys(int num_key1, unsigned char *k1,
									Index<L2<unsigned char> > index,
									double ratio, int max_pts_visit)
{
	std::vector<KeypointMatch> matches;
	SearchParams searchParams(max_pts_visit);
	for (size_t i = 0; i < num_key1; i++)
	{
		Matrix<unsigned char> query(k1 + 128 * i, 1, 128);

		Matrix<int> indices(new int[2], 1, 2);
		Matrix<float> dists(new float[2], 1, 2);

		index.knnSearch(query, indices, dists, 2, searchParams);
		if (((double)dists[0][0]) < ratio * ratio * ((double)dists[0][1]))
		{	
			matches.push_back(KeypointMatch(i, indices[0][0]));
		}

	}
		
	return matches;
}

