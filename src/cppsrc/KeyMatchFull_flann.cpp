#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <ctype.h>
#include <vector>
#include <algorithm>

#include "keys2a_flann.h"

int ReadFileList(char* list_in, std::vector<std::string>& key_files) {
	FILE* fp;

	if ((fp = fopen(list_in, "r")) == NULL) {
		printf("Error opening file %s for reading.\n", list_in);
		return 1;
	}

	char buf[512], *start;
	while (fgets(buf, 512, fp)) {
		// Remove trailing new-line
		if (buf[strlen(buf) - 1] == '\n') buf[strlen(buf) - 1] = '\0';

		// Find first non-space character
		start = buf;
		while (isspace(*start)) start++;

		// Skip empty lines
		if (strlen(start) == 0) continue;

		// Append file-name to key_files
		key_files.push_back(std::string(buf));
	}

	// Check we found input files
	if (key_files.size() == 0) {
		printf("No input files found in %s.\n", list_in);
		return 1;
	}

	return 0;
}

int main(int argc, char* argv[])
{
	char *list_in;
	char *file_out;
	double ratio;

	double ccd_width = 6.17;
	double f_mm = 3.6;
	
	if (argc != 3 && argc != 4)
	{
		printf("Usage: %s <list.txt> <outfile> [window_radius]\n", argv[0]);
		return EXIT_FAILURE;
	}

	list_in = argv[1];
	ratio = 0.6;
	file_out = argv[2];

	int window_radius = -1;
	if (argc == 4)
	{
		window_radius = atoi(argv[3]);
	}

	clock_t start = clock();

	/* Read the list of files*/
	std::vector<std::string> key_files;
	if (ReadFileList(list_in, key_files) != 0)	return EXIT_FAILURE;

	FILE *f;
	if ((f = fopen(file_out, "w")) == NULL)
	{
		printf("Could not open %s for writing.\n", file_out);
		return EXIT_FAILURE;
	}

	int num_images = (int)key_files.size();

	std::vector<unsigned char*> keys(num_images);
	std::vector<int> num_keys(num_images);

	/* Read all keys */
	for (size_t i = 0; i < num_images; i++)
	{
		keys[i] = NULL;
		num_keys[i] =  ReadKeyFile(key_files[i].c_str(), &keys[i]);
	}

	clock_t end = clock();
	printf("[KeyMatchFull] Reading keys took %0.3fs\n",
		(end - start) / ((double)CLOCKS_PER_SEC));
	
	
	for (size_t i = 0; i < num_images; i++)
	{
		if (num_keys[i] ==0)
			continue;

		printf("[KeyMatchFull] Matching to image %d\n", i);
		
		start = clock();

		/* Create a kdindex from the keys */
		Matrix<unsigned char> features(keys[i], num_keys[i], 128);
		Index<L2<unsigned char> > index(features, KDTreeIndexParams(1));
		index.buildIndex();


		/* Compute the start index */
		int start_idx = 0;
		if (window_radius > 0)
			start_idx = std::max((int)(i - window_radius), 0);

		for (size_t j = start_idx; j < i; j++)
		{
	
			if (num_keys[j] == 0)	continue;

			/* Compute likely matches between two sets of ketpoints */
		//	std::vector<KeypointMatch> matches = MatchKeys(num_keys[j], keys[j], index, ratio);


			std::vector<KeypointMatch> matches;
			SearchParams searchParams(10);
			for(size_t k = 0; k < num_keys[j]; k++)
			{
				Matrix<unsigned char> query(keys[j] + 128*k,1,128);
				std::vector<std::vector<int> > indices;
				std::vector<std::vector<float> > dists;
				index.knnSearch(query, indices, dists, 2, searchParams);
				if(((double) dists[0][0]) < ratio * ratio * ((double) dists[0][1]))
					matches.push_back(KeypointMatch(k,indices[0][0]));
			}
			int num_matches = (int)matches.size();

			if (num_matches >= 16)
			{
				/* Write the pair */
				fprintf(f, "%d %d\n", j, i);

				/* Write the number of matches */
				fprintf(f, "%d\n", (int)matches.size());

				for (size_t k = 0; k < num_matches; k++)
				{
					fprintf(f, "%d %d\n",
						matches[k].m_idx1, matches[k].m_idx2);
				}
			}
		}
		end = clock();
		printf("[KeyMatchFull] Matching took %0.3fs\n",
			(end - start) / (double)CLOCKS_PER_SEC);
		fflush(stdout);

	}

	/* Free keypoints */
	for (size_t i = 0; i < num_images; i++)
	{
		if (keys[i] != NULL)
			delete[] keys[i];
	}

	fclose(f);
	return EXIT_SUCCESS;
}

