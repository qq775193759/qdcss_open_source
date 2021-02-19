#include "SuperVoxel_MatricIMSLIC.h"

#include <iostream>
#include <queue>
#include <algorithm>
#include <cmath>
#include <time.h>

using namespace cv;
using namespace std;


SuperVoxel_MatricIMSLIC::SuperVoxel_MatricIMSLIC()
{
}


SuperVoxel_MatricIMSLIC::~SuperVoxel_MatricIMSLIC()
{
}

void SuperVoxel_MatricIMSLIC::Matric_DoSupervoxelSegmentation(const int pK, const double pC_xy, const double pC_z, const int iters)
{
	K = pK;
	C_xy = pC_xy;
	C_z = pC_z;

	klabels = vector<int>(sz_xyz, -1);
	hash = vector<int>(sz_xyz, 0);

	//framework
	Matric_MakeMatric();
	Matric_kmeanspp();
	for(int i=1;i<iters;i++)
	{
		Matric_Voronoi();
		Matric_Center();
	}
	Matric_Voronoi();
	Matric_Connectivity();
}

void SuperVoxel_MatricIMSLIC::Matric_kmeanspp()
{
	distvec = vector<double>(sz_xyz, 1e9);
	kseedx.resize(K);
	kseedy.resize(K);
	kseedz.resize(K);

	for(int k=0;k<K;k++)
	{
		double sum = 0;
		for(int i=0;i<sz_xyz;i++)
			sum += distvec[i]*distvec[i];
		double rand_01 = rand()/(RAND_MAX+0.1);
		double rand_val = rand_01*sum;
		sum=0;
		int rand_rank = 0;
		for(int i=0;i<sz_xyz;i++)
		{
			sum += distvec[i]*distvec[i];
			if(sum > rand_val)
			{
				rand_rank = i;
				break;
			}
		}
		kseedx[k] = rand_rank%m_width;
		kseedy[k] = (rand_rank%sz_xy)/m_width;
		kseedz[k] = rand_rank/sz_xy;
		int y1 = max(0,			kseedy[k]-offsety);
		int y2 = min(m_height,	kseedy[k]+offsety);
		int x1 = max(0,			kseedx[k]-offsetx);
		int x2 = min(m_width,	kseedx[k]+offsetx);
		int z1 = max(0,			kseedz[k]-offsetz);
		int z2 = min(m_frame,	kseedz[k]+offsetz);
		Matric_Dynamic_Dist_Without_Segment(kseedx[k], kseedy[k], kseedz[k], x1, y1, z1, x2, y2, z2, 0);
	}
}

void SuperVoxel_MatricIMSLIC::Matric_Voronoi()
{
	distvec = vector<double>(sz_xyz, DBL_MAX);

	for( int k = 0; k < K; k++ )
	{
		int y1 = max(0,			kseedy[k]-offsety);
		int y2 = min(m_height,	kseedy[k]+offsety);
		int x1 = max(0,			kseedx[k]-offsetx);
		int x2 = min(m_width,	kseedx[k]+offsetx);
		int z1 = max(0,			kseedz[k]-offsetz);
		int z2 = min(m_frame,	kseedz[k]+offsetz);

		Matric_Dynamic_Dist_Without_Segment(kseedx[k], kseedy[k], kseedz[k], x1, y1, z1, x2, y2, z2, k);
	}
}

void SuperVoxel_MatricIMSLIC::Matric_Center()
{
	vector<double> clustersizex(K, 0);
	vector<double> clustersizey(K, 0);
	vector<double> clustersizez(K, 0);
	vector<double> sigmax(K, 0);
	vector<double> sigmay(K, 0);
	vector<double> sigmaz(K, 0);
	for(int h=0;h<m_frame;h++)
		for(int i=0;i<m_height;i++)
			for(int j=0;j<m_width;j++)
			{
				int p0 = h*sz_xy + i*m_width + j;
				int current_label = klabels[p0];
				if(current_label >=0)
				{
					double powerx = 1;//1+distvec[p0];
					double powery = 1;//1+distvec[p0];
					double powerz = 1;//1+distvec[p0];
					sigmax[current_label] += j * powerx;
					sigmay[current_label] += i * powery;
					sigmaz[current_label] += h * powerz;
					clustersizex[current_label] += powerx;
					clustersizey[current_label] += powery;
					clustersizez[current_label] += powerz;
				}
			}
	for( int k = 0; k < K; k++ )
	{
		kseedx[k] = int(sigmax[k]/clustersizex[k] + 0.5);
		kseedy[k] = int(sigmay[k]/clustersizex[k] + 0.5);
		kseedz[k] = int(sigmaz[k]/clustersizex[k] + 0.5);
	}
}

void SuperVoxel_MatricIMSLIC::Matric_Connectivity()
{
	vector<int> nlabels(sz_xyz, -1);
	vector<int> hash(sz_xyz, 0);
	for(int k=0;k<K;k++)
	{
		queue<Voxel> pq_pixel;
		pq_pixel.push(Voxel(kseedx[k], kseedy[k], kseedz[k]));
		while(!pq_pixel.empty())
		{
			Voxel tmp_pixel = pq_pixel.front();
			int tx = tmp_pixel.x;
			int ty = tmp_pixel.y;
			int tz = tmp_pixel.z;
			pq_pixel.pop();
			nlabels[tz*sz_xy+ty*m_width+tx] = k;
			hash[tz*sz_xy+ty*m_width+tx] = 1;
			//go through all neighbor
			for(int i=0;i<10;i++)
			{
				int txi = tx + nb_dx10[i];
				int tyi = ty + nb_dy10[i];
				int tzi = tz + nb_dz10[i];
				if(txi < 0 || txi >= m_width || tyi < 0 || tyi >= m_height || tzi < 0 || tzi >= m_frame)
					continue;
				if(klabels[tzi*sz_xy+tyi*m_width+txi] != k)
					continue;
				if(hash[tzi*sz_xy+tyi*m_width+txi] == 1)
					continue;
				if(hash[tzi*sz_xy+tyi*m_width+txi] == 0)
				{
					pq_pixel.push(Voxel(txi, tyi, tzi));
					hash[tzi*sz_xy+tyi*m_width+txi] = 2;
				}
			}
		}
	}
	int co = 0;
	for(int i=0;i<sz_xyz;i++)
	{
		if(nlabels[i] == -1)
		{
			co++;
			int tmplabel;
			if((i%m_width)>0)
				tmplabel = nlabels[i-1];
			else if(((i%sz_xy)/m_width)>0)
				tmplabel = nlabels[i-m_width];
			else if((i/sz_xy)>0)
				tmplabel = nlabels[i-sz_xy];
			else
				tmplabel = 0;
			queue<Voxel> pq_pixel;
			pq_pixel.push(Voxel(i%m_width, (i%sz_xy)/m_width, i/sz_xy));
			while(!pq_pixel.empty())
			{
				Voxel tmp_pixel = pq_pixel.front();
				int tx = tmp_pixel.x;
				int ty = tmp_pixel.y;
				int tz = tmp_pixel.z;
				pq_pixel.pop();
				nlabels[tz*sz_xy+ty*m_width+tx] = tmplabel;
				hash[tz*sz_xy+ty*m_width+tx] = 1;
				//go through all neighbor
				for(int i=0;i<10;i++)
				{
					int txi = tx + nb_dx10[i];
					int tyi = ty + nb_dy10[i];
					int tzi = tz + nb_dz10[i];
					if(txi < 0 || txi >= m_width || tyi < 0 || tyi >= m_height || tzi < 0 || tzi >= m_frame)
						continue;
					if(nlabels[tzi*sz_xy+tyi*m_width+txi] != -1)
						continue;
					if(hash[tzi*sz_xy+tyi*m_width+txi] == 1)
						continue;
					if(hash[tzi*sz_xy+tyi*m_width+txi] == 0)
					{
						pq_pixel.push(Voxel(txi, tyi, tzi));
						hash[tzi*sz_xy+tyi*m_width+txi] = 2;
					}
				}
			}
		}
	}
}





double SuperVoxel_MatricIMSLIC::Matric_Dist(int x1, int y1, int z1, int x2, int y2, int z2)//for conputing dist
{
	if(x1 < 0 || x2 < 0 || y1 < 0 || y2 < 0 || x1 >= m_width || x2 >= m_width || y1 >= m_height || y2>= m_height)
		return 1e6;
	if(z1 < 0 || z2 < 0 || z1 >= m_frame || z2 >= m_frame)
		return 1e6;
	int p1 = z1*sz_xy + y1*m_width + x1;
	int p2 = z2*sz_xy + y2*m_width + x2;


	double l1 = m_lvec[p1];
	double a1 = m_avec[p1];
	double b1 = m_bvec[p1];
	double l2 = m_lvec[p2];
	double a2 = m_avec[p2];
	double b2 = m_bvec[p2];
	double dist = (l1 - l2)*(l1 - l2) + (a1 - a2)*(a1 - a2) + (b1 - b2)*(b1 - b2);

	//return (sqrt(dist) + sqrt(invwt));//dist = sqrt(dist) + sqrt(distxy*invwt);//this is more exact
	double dist_xy = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
	double dist_z = (z1-z2)*(z1-z2);
	return sqrt(dist + C_xy*dist_xy + C_z*dist_z);
}

void SuperVoxel_MatricIMSLIC::Matric_MakeMatric()
{
	//srand(19960211);
	//srand(19961204);
	//srand(20181116);
	srand((int)time(0));rand();
	//neighbor dx
	nb_dx10[0]=nb_dx10[2]=nb_dx10[8]=nb_dx10[9]=0;
	nb_dx10[4]=nb_dx10[1]=nb_dx10[5]=1;
	nb_dx10[7]=nb_dx10[3]=nb_dx10[6]=-1;
	nb_dx[0]=nb_dx[2]=nb_dx[8]=nb_dx[9]=nb_dx[11]=nb_dx[17]=nb_dx[18]=nb_dx[20]=0;
	nb_dx[4]=nb_dx[1]=nb_dx[5]=nb_dx[13]=nb_dx[10]=nb_dx[14]=nb_dx[22]=nb_dx[19]=nb_dx[23]=1;
	nb_dx[7]=nb_dx[3]=nb_dx[6]=nb_dx[16]=nb_dx[12]=nb_dx[15]=nb_dx[25]=nb_dx[21]=nb_dx[24]=-1;
	//neighbor dy
	nb_dx10[3]=nb_dx10[1]=nb_dx10[8]=nb_dx10[9]=0;
	nb_dx10[7]=nb_dx10[0]=nb_dx10[4]=1;
	nb_dx10[6]=nb_dx10[2]=nb_dx10[5]=-1;
	nb_dy[3]=nb_dy[1]=nb_dy[12]=nb_dy[8]=nb_dy[10]=nb_dy[21]=nb_dy[17]=nb_dy[19]=0;
	nb_dy[7]=nb_dy[0]=nb_dy[4]=nb_dy[16]=nb_dy[9]=nb_dy[13]=nb_dy[25]=nb_dy[18]=nb_dy[22]=1;
	nb_dy[6]=nb_dy[2]=nb_dy[5]=nb_dy[15]=nb_dy[11]=nb_dy[14]=nb_dy[24]=nb_dy[20]=nb_dy[23]=-1;
	//neighbor dz
	nb_dx10[0]=nb_dx10[1]=nb_dx10[2]=nb_dx10[3]=nb_dx10[4]=nb_dx10[5]=nb_dx10[6]=nb_dx10[7]=0;
	nb_dx10[8]=1;
	nb_dx10[9]=-1;
	nb_dz[0]=nb_dz[1]=nb_dz[2]=nb_dz[3]=nb_dz[4]=nb_dz[5]=nb_dz[6]=nb_dz[7]=0;
	nb_dz[8]=nb_dz[9]=nb_dz[10]=nb_dz[11]=nb_dz[12]=nb_dz[13]=nb_dz[14]=nb_dz[15]=nb_dz[16]=1;
	nb_dz[17]=nb_dz[18]=nb_dz[19]=nb_dz[20]=nb_dz[21]=nb_dz[22]=nb_dz[23]=nb_dz[24]=nb_dz[25]=-1;
	//offset
	int tK = sqrt(K/7);
	offsetx = m_width/tK;
	offsety = m_height/tK;
	offsetz = 20;
	/*for(int h=0;h<m_frame;h++)
		for(int i=0;i<m_height;i++)
			for(int j=0;j<m_width;j++)//width == x
				for(int k=0;k<10;k++)
				{
					double tmp_dist = Matric_Dist(j, i, h, j+nb_dx[k], i+nb_dy[k], h+nb_dz[k]);
					distNeighbor[k].push_back(tmp_dist);
				}*/
}


void SuperVoxel_MatricIMSLIC::Matric_Dynamic_Dist_Without_Segment(
	int x0, int y0, int z0, int x1, int y1, int z1, int x2, int y2, int z2, int label_no)
{
	for(int h=z1;h<z2;h++)
		for(int i=y1;i<y2;i++)
			for(int j=x1;j<x2;j++)
			{
				hash[h*sz_xy+i*m_width+j] = 0;
			}
	distvec[z0*sz_xy + y0*m_width + x0] = 0;
	queue<Voxel> pq_voxel;
	pq_voxel.push(Voxel(x0, y0, z0));
	while(!pq_voxel.empty())
	{
		Voxel tmp_voxel = pq_voxel.front();
		pq_voxel.pop();
		int tx = tmp_voxel.x;
		int ty = tmp_voxel.y;
		int tz = tmp_voxel.z;
		double tdist = distvec[tz*sz_xy+ty*m_width+tx];
		hash[tz*sz_xy+ty*m_width+tx] = 1;
		klabels[tz*sz_xy+ty*m_width+tx] = label_no;
		//go through all neighbor
		for(int i=0;i<26;i++)
		{
			int txi = tx + nb_dx[i];
			int tyi = ty + nb_dy[i];
			int tzi = tz + nb_dz[i];
			if(txi < x1 || txi >= x2 || tyi < y1 || tyi >= y2 || tzi < z1 || tzi >= z2)
				continue;
			if(hash[tzi*sz_xy+tyi*m_width+txi] == 1)
				continue;
			//double tmp_dist = distNeighbor[i][tz*sz_xy+ty*m_width+tx];
			double tmp_dist = Matric_Dist(tx, ty, tz, txi, tyi, tzi);
			if((tdist + tmp_dist) < distvec[tzi*sz_xy+tyi*m_width+txi])
			{
				distvec[tzi*sz_xy+tyi*m_width+txi] = tdist + tmp_dist;
				if(hash[tzi*sz_xy+tyi*m_width+txi] == 0)
				{
					pq_voxel.push(Voxel(txi, tyi, tzi));
					hash[tzi*sz_xy+tyi*m_width+txi] = 2;
				}
			}
		}
	}
}


void SuperVoxel_MatricIMSLIC::DrawContoursAroundSegments(int frame_i, const string out_name)
{
	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};
	vector<bool> istaken(sz_xy, false);
	int mainindex(0);
	for( int j = 0; j < m_height; j++ )
	{
		for( int k = 0; k < m_width; k++ )
		{
			int np(0);
			for( int i = 0; i < 8; i++ )
			{
				int x = k + dx8[i];
				int y = j + dy8[i];

				if( (x >= 0 && x < m_width) && (y >= 0 && y < m_height) )
				{
					int index = y*m_width + x;
					{
						if( istaken[index] == false)
							if( klabels[frame_i*sz_xy + mainindex] != klabels[frame_i*sz_xy + index] ) np++;
					}
				}
			}
			//if( np > 1 )//change to 2 or 3 for thinner lines
			if( np > 0 )
			{
				mat_vec[frame_i].at<Vec3b>(j, k)[0] = 0;
				mat_vec[frame_i].at<Vec3b>(j, k)[1] = 128;
				mat_vec[frame_i].at<Vec3b>(j, k)[2] = 0;
				istaken[mainindex] = true;
			}
			mainindex++;
		}
	}
	for(int k=0;k<K;k++)
	{
		if(kseedz[k] == frame_i)
		{
			mat_vec[frame_i].at<Vec3b>(kseedy[k], kseedx[k])[0] = 0;
			mat_vec[frame_i].at<Vec3b>(kseedy[k], kseedx[k])[1] = 0;
			mat_vec[frame_i].at<Vec3b>(kseedy[k], kseedx[k])[2] = 128;
		}
	}
	imwrite(out_name, mat_vec[frame_i]);
}


void SuperVoxel_MatricIMSLIC::save_line_res(const string out_path)
{
	for(int i=0;i<m_frame;i++)
	{
		string out_name = out_path + "\\" + to_string(i) + ".bmp";
		DrawContoursAroundSegments(i, out_name);
	}
}

void SuperVoxel_MatricIMSLIC::save_color_res(const string out_path)
{
	vector<uchar> rcolor[3];
	for(int i=0;i<3;i++)
		for(int k=0;k<K;k++)
		{
			uchar tmp_color = rand()%256;
			rcolor[i].push_back(tmp_color);
		}
	for(int h=0;h<m_frame;h++)
	{
		string out_name = out_path + "\\" + pngname[h];
		for(int i=0;i<m_height;i++)
			for(int j=0;j<m_width;j++)
				for(int k=0;k<3;k++)
				{
					mat_vec[h].at<Vec3b>(i, j)[k] = rcolor[k][klabels[h*sz_xy + i*m_width + j]];
				}
		imwrite(out_name, mat_vec[h]);
	}
}


