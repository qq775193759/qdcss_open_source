#include "MainProcess.h"
#include <time.h>
#include <io.h>
#include <thread>
using namespace std;
using namespace cv;

void label2csv(const string out_name, int* labels, int width, int height)
{
	ofstream fout(out_name);
	int sz = width * height;
	for(int i=0;i<sz;i++)
	{
		fout<<labels[i]<<',';
		//fout<<labels[i]<<' ';
		if((i%width)==(width-1))
			fout<<endl;
	}
	fout.close();
}

string myStrReplace(const string s, const string s1, const string s2)
{
	string res = s;
	int pos = res.find(s1);
	if(pos >= res.size())
		return res;
	res.replace(pos, s1.size(), s2);
	return res;
}

void getFiles( string path, vector<string>& files )  
{   
#if _WIN64   
	__int64     hFile   =   0;  
#else
    long   hFile   =   0;  
#endif
    struct _finddata_t fileinfo;  
    string p;  
    if((hFile = _findfirst(p.assign(path).append("\\*").c_str(),&fileinfo)) !=  -1)  
    {  
        do  
        {   
			if(!(fileinfo.attrib &  _A_SUBDIR))
                files.push_back(p.assign(path).append("\\").append(fileinfo.name) );  
        }
		while(_findnext(hFile, &fileinfo)  == 0);  
        _findclose(hFile);  
    }  
}

void RGB2XYZ(const int&	sR, const int&	sG, const int&	sB, double&	X, double& Y, double& Z)
{
	double R = sR/255.0;
	double G = sG/255.0;
	double B = sB/255.0;

	double r, g, b;

	if(R <= 0.04045)	r = R/12.92;
	else				r = pow((R+0.055)/1.055,2.4);
	if(G <= 0.04045)	g = G/12.92;
	else				g = pow((G+0.055)/1.055,2.4);
	if(B <= 0.04045)	b = B/12.92;
	else				b = pow((B+0.055)/1.055,2.4);

	X = r*0.4124564 + g*0.3575761 + b*0.1804375;
	Y = r*0.2126729 + g*0.7151522 + b*0.0721750;
	Z = r*0.0193339 + g*0.1191920 + b*0.9503041;
}

void RGB2LAB(const int& sR, const int& sG, const int& sB, double& lval, double& aval, double& bval)
{
	//------------------------
	// sRGB to XYZ conversion
	//------------------------
	double X, Y, Z;
	RGB2XYZ(sR, sG, sB, X, Y, Z);

	//------------------------
	// XYZ to LAB conversion
	//------------------------
	double epsilon = 0.008856;	//actual CIE standard
	double kappa   = 903.3;		//actual CIE standard

	double Xr = 0.950456;	//reference white
	double Yr = 1.0;		//reference white
	double Zr = 1.088754;	//reference white

	double xr = X/Xr;
	double yr = Y/Yr;
	double zr = Z/Zr;

	double fx, fy, fz;
	if(xr > epsilon)	fx = pow(xr, 1.0/3.0);
	else				fx = (kappa*xr + 16.0)/116.0;
	if(yr > epsilon)	fy = pow(yr, 1.0/3.0);
	else				fy = (kappa*yr + 16.0)/116.0;
	if(zr > epsilon)	fz = pow(zr, 1.0/3.0);
	else				fz = (kappa*zr + 16.0)/116.0;

	lval = 116.0*fy-16.0;
	aval = 500.0*(fx-fy);
	bval = 200.0*(fy-fz);
}

void superPixel_rgb_iter_once(const string pic_name, const string out_name, int K, int iter)
{
	//read picture
	Mat mat = imread(pic_name);
	int width = mat.cols;
	int height = mat.rows;
	int sz = width * height;
	int* labels = new int[sz];
	unsigned int* img = new unsigned int[sz];
	for(int i=0;i<height;i++)
		for(int j=0;j<width;j++)//width == x
		{
			img[i*width+j] = mat.at<Vec3b>(i, j)[2];
			img[i*width+j] = (img[i*width+j]<<8) + mat.at<Vec3b>(i, j)[1];
			img[i*width+j] = (img[i*width+j]<<8) + mat.at<Vec3b>(i, j)[0];
		}

	MSLIC MSLIC;
	MSLIC.ifdraw = 1;
	MSLIC.iteration = iter;

	MSLIC.Matric_DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(img, width, height, labels, K, 20);

	MSLIC.DrawContoursAroundSegments(img, labels, width, height, 0);
	for(int i=0;i<height;i++)
		for(int j=0;j<width;j++)//width == x
		{
			mat.at<Vec3b>(i, j)[2] = (img[i*width+j] >> 16) & 0xFF;
			mat.at<Vec3b>(i, j)[1] = (img[i*width+j] >> 8) & 0xFF;
			mat.at<Vec3b>(i, j)[0] = img[i*width+j] & 0xFF;
		}
	imwrite(out_name, mat);

	mat.release();
}

void superPixel_rgb_once(const string pic_name, const string out_name, int m_compactness, int K)
{
	//read picture
	Mat mat = imread(pic_name);
	int width = mat.cols;
	int height = mat.rows;
	int sz = width * height;
	int* labels = new int[sz];
	unsigned int* img = new unsigned int[sz];
	for(int i=0;i<height;i++)
		for(int j=0;j<width;j++)//width == x
		{
			img[i*width+j] = mat.at<Vec3b>(i, j)[2];
			img[i*width+j] = (img[i*width+j]<<8) + mat.at<Vec3b>(i, j)[1];
			img[i*width+j] = (img[i*width+j]<<8) + mat.at<Vec3b>(i, j)[0];
		}

	MSLIC MSLIC;
	MSLIC.ifdraw = 1;
	MSLIC.iteration = 10;

	double dur;
	clock_t start,end;
	start = clock();
	MSLIC.Matric_DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(img, width, height, labels, K, m_compactness);
	end = clock();
	dur = (double)(end - start);
	printf("Use Time:%f\n",(dur/CLOCKS_PER_SEC));

	cout<<out_name<<endl;
	MSLIC.DrawContoursAroundSegments(img, labels, width, height, 0);
	for(int i=0;i<height;i++)
		for(int j=0;j<width;j++)//width == x
		{
			mat.at<Vec3b>(i, j)[2] = (img[i*width+j] >> 16) & 0xFF;
			mat.at<Vec3b>(i, j)[1] = (img[i*width+j] >> 8) & 0xFF;
			mat.at<Vec3b>(i, j)[0] = img[i*width+j] & 0xFF;
		}
	//imwrite(out_name, mat);
	cout<<myStrReplace(out_name, ".bmp", ".csv")<<endl;
	label2csv(myStrReplace(out_name, ".bmp", ".csv"), labels, width, height);
	//cout<<myStrReplace(out_name, ".bmp", ".txt")<<endl;
	//label2csv(myStrReplace(out_name, ".bmp", ".txt"), labels, width, height);

	mat.release();
}

void superPixel_rgbd_autodp_once(const string pic_name, const string out_name, double depth_power, int m_compactness, int K)
{
	//read picture
	Mat mat = imread(pic_name);
	int width = mat.cols;
	int height = mat.rows;
	int sz = width * height;
	int* labels = new int[sz];
	unsigned int* img = new unsigned int[sz];
	for(int i=0;i<height;i++)
		for(int j=0;j<width;j++)//width == x
		{
			img[i*width+j] = mat.at<Vec3b>(i, j)[2];
			img[i*width+j] = (img[i*width+j]<<8) + mat.at<Vec3b>(i, j)[1];
			img[i*width+j] = (img[i*width+j]<<8) + mat.at<Vec3b>(i, j)[0];
		}

	MSLIC MSLIC;
	MSLIC.ifdraw = 1;
	MSLIC.iteration = 10;
	MSLIC.depth_power = depth_power;
	MSLIC.depth_flag = 1;

	//rgbd 
	string depth_name = pic_name;
	depth_name = myStrReplace(depth_name, "images", "depth");
	//depth_name = myStrReplace(depth_name, "png", "txt");
	cout<<depth_name<<endl;
	Mat depth_mat = imread(depth_name);
	for(int i=0;i<height;i++)
		for(int j=0;j<width;j++)//width == x
		{
			MSLIC.depth_vec.push_back(depth_mat.at<unsigned short>(i, j));
		}

	double dur;
	clock_t start,end;
	start = clock();
	MSLIC.Matric_DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(img, width, height, labels, K, m_compactness);
	end = clock();
	dur = (double)(end - start);
	printf("Use Time:%f\n",(dur/CLOCKS_PER_SEC));

	cout<<out_name<<endl;
	MSLIC.DrawContoursAroundSegments(img, labels, width, height, 0);
	for(int i=0;i<height;i++)
		for(int j=0;j<width;j++)//width == x
		{
			mat.at<Vec3b>(i, j)[2] = (img[i*width+j] >> 16) & 0xFF;
			mat.at<Vec3b>(i, j)[1] = (img[i*width+j] >> 8) & 0xFF;
			mat.at<Vec3b>(i, j)[0] = img[i*width+j] & 0xFF;
		}
	imwrite(out_name, mat);


	cout<<myStrReplace(out_name, ".bmp", ".csv")<<endl;
	label2csv(myStrReplace(out_name, ".bmp", ".csv"), labels, width, height);
	mat.release();
}

void superPixel_rgbd_iter_once(const string pic_name, const string depth_name, const string out_name, int K, int iter)
{
	//read picture
	Mat mat = imread(pic_name);
	int width = mat.cols;
	int height = mat.rows;
	int sz = width * height;
	int* labels = new int[sz];
	unsigned int* img = new unsigned int[sz];
	for(int i=0;i<height;i++)
		for(int j=0;j<width;j++)//width == x
		{
			img[i*width+j] = mat.at<Vec3b>(i, j)[2];
			img[i*width+j] = (img[i*width+j]<<8) + mat.at<Vec3b>(i, j)[1];
			img[i*width+j] = (img[i*width+j]<<8) + mat.at<Vec3b>(i, j)[0];
		}

	MSLIC MSLIC;
	MSLIC.ifdraw = 1;
	MSLIC.iteration = iter;
	MSLIC.depth_power = 1e-4;
	MSLIC.depth_flag = 1;

	//rgbd 
	Mat depth_mat = imread(depth_name);
	for(int i=0;i<height;i++)
		for(int j=0;j<width;j++)//width == x
		{
			MSLIC.depth_vec.push_back(depth_mat.at<unsigned short>(i, j));
		}

	MSLIC.Matric_DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(img, width, height, labels, K, 20);

	MSLIC.DrawContoursAroundSegments(img, labels, width, height, 0);
	for(int i=0;i<height;i++)
		for(int j=0;j<width;j++)//width == x
		{
			mat.at<Vec3b>(i, j)[2] = (img[i*width+j] >> 16) & 0xFF;
			mat.at<Vec3b>(i, j)[1] = (img[i*width+j] >> 8) & 0xFF;
			mat.at<Vec3b>(i, j)[0] = img[i*width+j] & 0xFF;
		}
	imwrite(out_name, mat);
}


void superPixel_rgb(const string pic_path, const string out_path, int m_compactness, int K)
{
	vector<string> files;
	getFiles(pic_path, files);
	for(int i=0;i<files.size();i++)
	{
		string pic_name = files[i];
		string out_name = myStrReplace(pic_name, pic_path, out_path);
		out_name = myStrReplace(out_name, ".png", ".bmp");
		out_name = myStrReplace(out_name, ".jpg", ".bmp");
		cout<<pic_name<<endl;
		superPixel_rgb_once(pic_name, out_name, m_compactness, K);
	}
}

void superPixel_rgb_iter(const string pic_path, const string out_path, int K, int iter)
{
	vector<string> files;
	getFiles(pic_path, files);
	for(int i=0;i<files.size();i++)
	{
		string pic_name = files[i];
		string out_name = myStrReplace(pic_name, pic_path, out_path);
		out_name = myStrReplace(out_name, ".png", ".bmp");
		out_name = myStrReplace(out_name, ".jpg", ".bmp");
		cout<<pic_name<<endl;
		superPixel_rgb_iter_once(pic_name, out_name, K, iter);
	}
}

const int THREAD_MAX = 8;

void superPixel_rgb_iter_mt_func(const string pic_path, const string out_path, int K, int iter, int thread_i)
{
	vector<string> files;
	getFiles(pic_path, files);
	for(int i=thread_i;i<files.size();i+=THREAD_MAX)
	{
		string pic_name = files[i];
		string out_name = myStrReplace(pic_name, pic_path, out_path);
		out_name = myStrReplace(out_name, ".png", ".bmp");
		out_name = myStrReplace(out_name, ".jpg", ".bmp");
		superPixel_rgb_iter_once(pic_name, out_name, K, iter);
	}
}

void superPixel_rgb_iter_mt(const string pic_path, const string out_path, int K, int iter)
{
	thread t[THREAD_MAX];

	for(int i=0;i<THREAD_MAX;i++)
	{
		t[i] = thread(superPixel_rgb_iter_mt_func, pic_path, out_path, K, iter, i);
	}
	for(int i=0;i<THREAD_MAX;i++)
	{
		t[i].join();
	}
}

void superPixel_rgbd_autodp(const string pic_path, const string out_path, double depth_power, int m_compactness, int K)
{
	vector<string> files;
	getFiles(pic_path, files);
	for(int i=0;i<files.size();i++)
	{
		string pic_name = files[i];
		string out_name = myStrReplace(pic_name, pic_path, out_path);
		out_name = myStrReplace(out_name, ".png", ".bmp");
		out_name = myStrReplace(out_name, ".jpg", ".bmp");
		cout<<pic_name<<endl;
		superPixel_rgbd_autodp_once(pic_name, out_name, depth_power, m_compactness, K);
	}
}

void superPixel_rgbd_iter(const string pic_path, const string depth_path, const string out_path, int K, int iter)
{
	vector<string> files;
	getFiles(pic_path, files);
	for(int i=0;i<files.size();i++)
	{
		string pic_name = files[i];
		string depth_name = myStrReplace(pic_name, pic_path, depth_path);
		string out_name = myStrReplace(pic_name, pic_path, out_path);
		out_name = myStrReplace(out_name, ".png", ".bmp");
		out_name = myStrReplace(out_name, ".jpg", ".bmp");
		superPixel_rgbd_iter_once(pic_name, depth_name, out_name, K, iter);
	}
}

struct string_depth
{
	string pic;
	string depth;
	string out;
};

void superPixel_rgbd_iter_mt_func(const string_depth sd, int K, int iter, int thread_i)
{
	string pic_path = sd.pic;
	string depth_path = sd.depth;
	string out_path = sd.out;
	vector<string> files;
	getFiles(pic_path, files);
	for(int i=thread_i;i<files.size();i+=THREAD_MAX)
	{
		string pic_name = files[i];
		string depth_name = myStrReplace(pic_name, pic_path, depth_path);
		string out_name = myStrReplace(pic_name, pic_path, out_path);
		out_name = myStrReplace(out_name, ".png", ".bmp");
		out_name = myStrReplace(out_name, ".jpg", ".bmp");
		superPixel_rgbd_iter_once(pic_name, depth_name, out_name, K, iter);
	}
}

void superPixel_rgbd_iter_mt(const string pic_path, const string depth_path, const string out_path, int K, int iter)
{
	thread t[THREAD_MAX];

	for(int i=0;i<THREAD_MAX;i++)
	{
		string_depth tmp_sd;
		tmp_sd.pic = pic_path;
		tmp_sd.depth = depth_path;
		tmp_sd.out = out_path;
		t[i] = thread(superPixel_rgbd_iter_mt_func, tmp_sd, K, iter, i);
	}
	for(int i=0;i<THREAD_MAX;i++)
	{
		t[i].join();
	}
}

void superVoxel_readPath(const string pic_path, SuperVoxel_MatricIMSLIC& imslic)
{
	vector<string> files;
	getFiles(pic_path, files);
	cout<<"Frames: "<<files.size()<<endl;
	if(files.size() == 0)
		return;
	Mat tmp_mat = imread(files[0]);
	imslic.m_width = tmp_mat.cols;
	imslic.m_height = tmp_mat.rows;
	imslic.m_frame = files.size();
	imslic.sz_xy = imslic.m_width * imslic.m_height;
	imslic.sz_xyz = imslic.sz_xy * imslic.m_frame;
	for(int h=0;h<files.size();h++)
	{
		Mat mat = imread(files[h]);
		imslic.pngname.push_back(files[h].substr(files[h].find_last_of("\\") + 1));
		imslic.mat_vec.push_back(mat);
		//cout<<files[h]<<endl;
		//cout<<mat.size()<<endl;
		for(int i=0;i<imslic.m_height;i++)
			for(int j=0;j<imslic.m_width;j++)//width == x
			{
				double tl, ta, tb;
				RGB2LAB(mat.at<Vec3b>(i, j)[2], mat.at<Vec3b>(i, j)[1], mat.at<Vec3b>(i, j)[0], tl, ta, tb);
				imslic.m_lvec.push_back(tl);
				imslic.m_avec.push_back(ta);
				imslic.m_bvec.push_back(tb);
			}
	}
}

void superVoxel_rgb(const string pic_path, const string out_path, const int pK, const double pC_xy, const double pC_z, const int iters)
{
	SuperVoxel_MatricIMSLIC imslic;
	superVoxel_readPath(pic_path, imslic);

	double dur;
	clock_t start,end;
	start = clock();
	imslic.Matric_DoSupervoxelSegmentation(pK, pC_xy, pC_z, iters);
	end = clock();
	dur = (double)(end - start);
	printf("Use Time:%f\n",(dur/CLOCKS_PER_SEC));

	imslic.save_color_res(out_path);
}

void superVoxel_submit(const string pic_path, const string out_path, const int pK, const int iters)
{
	SuperVoxel_MatricIMSLIC imslic;
	superVoxel_readPath(pic_path, imslic);

	imslic.Matric_DoSupervoxelSegmentation(pK, 0.5, 0.5, iters);
	imslic.save_color_res(out_path);
}