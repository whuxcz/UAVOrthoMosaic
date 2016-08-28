#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <cstring>
#include <stdlib.h>

int main(int argc, char* argv[])
{
	if(argc < 4)
	{
		std::cout << "[Usage] 1.ply 2.ply ... n.ply final.ply" << std::endl;
		return -1;
	}	
	std::vector<char *> ply_file;
	char *ply_out = argv[argc-1];
	for(size_t i = 1; i < argc; i++)
	{
		ply_file.push_back(argv[i]);
	}	
	int pts_sum = 0;
	int ply_num = argc-2;
	std::ifstream iPly;
	std::ofstream oPly;
	std::string line;
	std::vector<int> pts_num;
	for (size_t i = 0; i < ply_num; i++)
	{
		iPly.open(ply_file[i]);
		if (!iPly)
		{
			std::cout << "wrong" << std::endl;
		}
		std::getline(iPly, line);
		if (line.find_first_of('p') != 0)
			return -1;

		std::getline(iPly, line);
		std::getline(iPly, line);

		char *p = (char *)line.data();
		strtok(p, " ");
		strtok(NULL, " ");

		int num = atoi((const char *)strtok(NULL, " "));
		pts_sum += num;
		pts_num.push_back(num);
		iPly.close();
	}
	std::cout << pts_sum << std::endl;
	oPly.open(ply_out);
	oPly << "ply\nformat ascii 1.0\nelement vertex " << pts_sum << "\nproperty float x\nproperty float y\nproperty float z\n"
		<< "property float nx\nproperty float ny\nproperty float nz\n" << "property uchar diffuse_red\nproperty uchar diffuse_green\nproperty uchar diffuse_blue\nend_header\n";
	for (size_t i = 0; i < ply_num; i++)
	{
		iPly.open(ply_file[i]);
		if (!iPly)
		{
			std::cout << "wrong" << std::endl;
		}
		for (size_t j = 0; j < 13; j++)
			std::getline(iPly, line);

		if (line.find_first_of('e') != 0)
		{
			std::cout << "wrong ply header" << std::endl;
			return -1;
		}

		for (size_t j = 0; j < pts_num[i]; j++)
		{
			std::getline(iPly, line);
			oPly << line << std::endl;;
		}

		iPly.close();
	}
	oPly.close();

	return 0;
}

