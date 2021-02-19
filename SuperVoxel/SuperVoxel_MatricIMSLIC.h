#pragma once
#include<vector>
#include <opencv2\opencv.hpp>
using namespace cv;
using namespace std;


class SuperVoxel_MatricIMSLIC
{
	vector<int> kseedx;
	vector<int> kseedy;
	vector<int> kseedz;
	vector<int> klabels;
	int K;
	double C_xy, C_z;

	//dist
	//vector<double> distNeighbor[10];
	int nb_dx[26];
	int nb_dy[26];
	int nb_dz[26];
	//  7 0 4     16  9 13    25 18 22
	//  3   1     12  8 10    21 17 19
	//  6 2 5     15 11 14    24 20 23
	int nb_dx10[10];
	int nb_dy10[10];
	int nb_dz10[10];
	//path length
	vector<double> distvec;
	vector<int> hash;

	//private tools
	double Matric_Dist(int x1, int y1, int z1, int x2, int y2, int z2);//for conputing dist
	void Matric_MakeMatric();
	void Matric_Dynamic_Dist_Without_Segment(
		int x0, int y0, int z0, int x1, int y1, int z1, int x2, int y2, int z2, int label_no);

	//main function
	int offsetx, offsety, offsetz;
	void Matric_kmeanspp();
	void Matric_Voronoi();
	void Matric_Center();
	void Matric_Connectivity();
public:
	SuperVoxel_MatricIMSLIC();
	~SuperVoxel_MatricIMSLIC();

	int m_width, m_height, m_frame;
	int sz_xy;
	int sz_xyz;
	vector<string> pngname;
	vector<double> m_lvec;
	vector<double> m_avec;
	vector<double> m_bvec;
	vector<Mat> mat_vec;


	void Matric_DoSupervoxelSegmentation(const int pK, const double pC_xy, const double pC_z, const int iters);


	void DrawContoursAroundSegments(int frame_i, const string out_name);
	void save_line_res(const string out_path);
	void save_color_res(const string out_path);
};



struct Voxel{
	int x;
	int y;
	int z;
	Voxel(int px, int py, int pz):x(px), y(py), z(pz){}
};