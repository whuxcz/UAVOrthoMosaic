#include <string>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>

int ModifyKey(const char *inKey, double mul, const char* outKey)
{
	std::ifstream iKey;
	std::ofstream oKey;
	iKey.open(inKey);
	oKey.open(outKey);
	std::string line;
	int key_num;
	int dim;
	iKey >> key_num >> dim;
	oKey << key_num << " " << dim << std::endl;
	std::getline(iKey, line);
	if (dim != 128)
		return -1;
	for (size_t i = 0; i < key_num; i++)
	{
		float x, y, scale, ori;
		iKey >> x >> y >> scale >> ori;
		std::getline(iKey, line);

		oKey << x * mul << " " << y * mul << " " << scale * mul << " " << ori << std::endl;
		for (size_t j = 0; j < 7; j++)
		{
			std::getline(iKey, line);
			oKey << line << std::endl;
		}
	}
	iKey.close();
	oKey.close();
	return 1;
}

int main(int argc, char* argv[])
{
	if (argc != 4)
	{
		printf("[Usage]: %s in_keyfile scale mod_keyfile\n", argv[0]);
		return -1;
	}

	char *inKey = argv[1];
	char *outKey = argv[3];
	double scale= atof(argv[2]);
	
	ModifyKey(inKey, scale, outKey);

	return 1;
}


