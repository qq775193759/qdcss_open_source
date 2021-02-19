#include "MainProcess.h"
#include <stdlib.h>
#include <time.h>

using namespace std;
using namespace cv;




int main(int argc, char **argv)
{
	double dur;
	clock_t start,end;
	start = clock();

	/*
	cout<<"...SuperPixels RGB Version..."<<endl;
	const int PARAMETER_NUM = 5;
	if(argc < PARAMETER_NUM)
	{
		cout<<"USAGE:  *.exe  input_fold  output_fold  number_of_superpixels  iter_max"<<endl;
		return 0;
	}
	superPixel_rgb_iter_mt(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]));
	*/

	/*
	cout<<"...SuperPixels RGBD Version..."<<endl;
	const int PARAMETER_NUM = 6;
	if(argc < PARAMETER_NUM)
	{
		cout<<"USAGE:  *.exe  rgb_image_fold  depth_image_fold  output_fold  number_of_superpixels  iter_max"<<endl;
		return 0;
	}
	superPixel_rgbg_iter_mt(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]));
	*/


	
	cout<<"...SuperVoxels Version..."<<endl;
	const int PARAMETER_NUM = 5;
	if(argc < PARAMETER_NUM)
	{
		cout<<"USAGE:  *.exe  input_fold  output_fold  number_of_supervoxels  iter_max"<<endl;
		return 0;
	}
	superVoxel_submit(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]));
	




	/*const int PARAMETER_NUM = 7;
	if(argc < PARAMETER_NUM)
		return 0;
	superVoxel_rgb(argv[1], argv[2], atoi(argv[3]), atof(argv[4]), atof(argv[5]), atoi(argv[6]));*/

	/*const int PARAMETER_NUM = 5;
	if(argc < PARAMETER_NUM)
		return 0;
	superPixel_rgb(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]));*/

	/*const int PARAMETER_NUM = 6;
	if(argc < PARAMETER_NUM)
		return 0;
	superPixel_rgbd_autodp(argv[1], argv[2], atof(argv[3]), atoi(argv[4]), atoi(argv[5]));*/

	end = clock();
	dur = (double)(end - start);
	printf("Total Time:  %f\n",(dur/CLOCKS_PER_SEC));

	cout<<"Done! You can check the results in the output folder now."<<endl;

	return 0;
}