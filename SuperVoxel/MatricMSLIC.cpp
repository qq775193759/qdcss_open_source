#include <cfloat>
#include <cmath>
#include <iostream>
#include <fstream>
#include "MSLIC.h"
#include <map>
#include <set>
#include <queue>
#include <time.h>

//====================================================================================================================================================
//======================================================================================================
//======================================================================================================
//======================================================================================================
//======================================================================================================
//======================================================================================================
//======================================================================================================
//======================================================================================================
//======================================================================================================
//======================================================================================================
// Matric Space method implement
// coded by YZP
// d(x, y) : straight stretch in parameter space with a anisotropy matric
//======================================================================================================


double MSLIC::Matric_Dist(int x1, int y1, int x2, int y2)
{
	if(x1 < 0 || x2 < 0 || y1 < 0 || y2 < 0 || x1 >= m_width || x2 >= m_width || y1 >= m_height || y2>= m_height)
		return 1e6;
	int p1 = y1*m_width + x1;
	int p2 = y2*m_width + x2;

	double l1 = m_lvec[p1];
	double a1 = m_avec[p1];
	double b1 = m_bvec[p1];
	double l2 = m_lvec[p2];
	double a2 = m_avec[p2];
	double b2 = m_bvec[p2];
	double dist = (l1 - l2)*(l1 - l2) + (a1 - a2)*(a1 - a2) + (b1 - b2)*(b1 - b2);

	//return (sqrt(dist) + sqrt(invwt));//dist = sqrt(dist) + sqrt(distxy*invwt);//this is more exact
	double dist_xy = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
	if(depth_flag == 0)
		return sqrt(dist + invwt*dist_xy);
	else
	{
		double dist_depth = depth_vec[p1]-depth_vec[p2];
		dist_depth = dist_depth * dist_depth;
		return sqrt(dist + invwt*dist_xy + depth_power*dist_depth);
	}
	//return 1;
}

void MSLIC::Matric_MakeMatric()
{
	//srand(19960211);
	//srand(19961204);
	srand(20181116);
	//srand((int)time(0));rand();
	//neighbor dx
	nb_dx[0]=nb_dx[2]=0;
	nb_dx[4]=nb_dx[1]=nb_dx[5]=1;
	nb_dx[7]=nb_dx[3]=nb_dx[6]=-1;
	//neighbor dy
	nb_dy[3]=nb_dy[1]=0;
	nb_dy[7]=nb_dy[0]=nb_dy[4]=1;
	nb_dy[6]=nb_dy[2]=nb_dy[5]=-1;
	// (i, j) -> (i+1, j) & (i, j) -> (i, j+1)
	//rgb dist -> matrix
	for(int i=0;i<m_height;i++)
		for(int j=0;j<m_width;j++)//width == x
			for(int k=0;k<8;k++)
			{
				double tmp_dist = Matric_Dist(j, i, j+nb_dx[k], i+nb_dy[k]);
				distNeighbor[k].push_back(tmp_dist);
			}
}

struct Pixel{
	int x;
	int y;
	Pixel(int px, int py):x(px), y(py){}
};

struct Pixel_dist{
	int x;
	int y;
	double dist;
	Pixel_dist(int px, int py, double pd):x(px), y(py), dist(pd)
	{
	}
	bool operator<(const Pixel_dist &tar)const
	{
		//minimum heap
		if(dist == tar.dist)
		{
			if(x == tar.x)
				return y < tar.y;
			return x < tar.x;
		}
		return dist < tar.dist;
	}
};

void MSLIC::Matric_Dijkstra_Dist(vector<double> &res, vector<int> &hash, int x0, int y0, int x1, int y1, int x2, int y2)
{
	for(int i=y1;i<y2;i++)
		for(int j=x1;j<x2;j++)
		{
			res[i*m_width+j] = DBL_MAX;
			hash[i*m_width+j] = 0;
		}
	res[y0*m_width + x0] = 0;
	set<Pixel_dist> pq_pixel;
	pq_pixel.insert(Pixel_dist(x0, y0, 0));
	while(!pq_pixel.empty())
	{
		Pixel_dist tmp_pixel = *pq_pixel.begin();
		pq_pixel.erase(pq_pixel.begin());
		int tx = tmp_pixel.x;
		int ty = tmp_pixel.y;
		double tdist = tmp_pixel.dist;
		hash[ty*m_width+tx] = 1;
		//go through all neighbor
		for(int i=0;i<8;i++)
		{
			int txi = tx + nb_dx[i];
			int tyi = ty + nb_dy[i];
			if(txi < x1 || txi >= x2 || tyi < y1 || tyi >= y2)
				continue;
			if(hash[tyi*m_width+txi] == 1)
				continue;
			if((tdist + distNeighbor[i][ty*m_width+tx]) < res[tyi*m_width+txi])
			{
				set<Pixel_dist>::iterator tmp_it = pq_pixel.find(Pixel_dist(txi, tyi, res[tyi*m_width+txi]));
				if(tmp_it != pq_pixel.end())
				{
					pq_pixel.erase(tmp_it);
				}
				res[tyi*m_width+txi] = tdist + distNeighbor[i][ty*m_width+tx];
				pq_pixel.insert(Pixel_dist(txi, tyi, res[tyi*m_width+txi]));
			}
		}
	}
}

//using dynamic programming to approximate dijstra
void MSLIC::Matric_Dynamic_Dist(vector<double> &res, vector<int> &hash, int x0, int y0, int x1, int y1, int x2, int y2)
{
	for(int i=y1;i<y2;i++)
		for(int j=x1;j<x2;j++)
		{
			res[i*m_width+j] = 1e9;
			hash[i*m_width+j] = 0;
		}
	res[y0*m_width + x0] = 0;
	queue<Pixel> pq_pixel;
	pq_pixel.push(Pixel(x0, y0));
	while(!pq_pixel.empty())
	{
		Pixel tmp_pixel = pq_pixel.front();
		pq_pixel.pop();
		int tx = tmp_pixel.x;
		int ty = tmp_pixel.y;
		double tdist = res[ty*m_width+tx];
		hash[ty*m_width+tx] = 1;
		//go through all neighbor
		for(int i=0;i<8;i++)
		{
			int txi = tx + nb_dx[i];
			int tyi = ty + nb_dy[i];
			if(txi < x1 || txi >= x2 || tyi < y1 || tyi >= y2)
				continue;
			if(hash[tyi*m_width+txi] == 1)
				continue;
			if(hash[tyi*m_width+txi] == 0)
			{
				pq_pixel.push(Pixel(txi, tyi));
				hash[tyi*m_width+txi] = 2;
			}
			res[tyi*m_width+txi] = min(tdist + distNeighbor[i][ty*m_width+tx], res[tyi*m_width+txi]);
			//res[tyi*m_width+txi] = min(tdist + Matric_Dist(tx, ty, txi, tyi), res[tyi*m_width+txi]);
		}
	}
}

//using dynamic programming to approximate dijstra
void MSLIC::Matric_Dynamic_Dist_Without_Segment(int*& klabels, vector<int> &hash, int x0, int y0, int x1, int y1, int x2, int y2, int label_no)
{
	for(int i=y1;i<y2;i++)
		for(int j=x1;j<x2;j++)
		{
			hash[i*m_width+j] = 0;
		}
	distvec[y0*m_width + x0] = 0;
	queue<Pixel> pq_pixel;
	pq_pixel.push(Pixel(x0, y0));
	while(!pq_pixel.empty())
	{
		Pixel tmp_pixel = pq_pixel.front();
		pq_pixel.pop();
		int tx = tmp_pixel.x;
		int ty = tmp_pixel.y;
		double tdist = distvec[ty*m_width+tx];
		hash[ty*m_width+tx] = 1;
		klabels[ty*m_width+tx] = label_no;
		//go through all neighbor
		for(int i=0;i<8;i++)
		{
			int txi = tx + nb_dx[i];
			int tyi = ty + nb_dy[i];
			if(txi < x1 || txi >= x2 || tyi < y1 || tyi >= y2)
				continue;
			if(hash[tyi*m_width+txi] == 1)
				continue;
			double tmp_dist = distNeighbor[i][ty*m_width+tx];
			//double tmp_dist = Matric_Dist(tx, ty, txi, tyi);
			if((tdist + tmp_dist) < distvec[tyi*m_width+txi])
			{
				distvec[tyi*m_width+txi] = tdist + tmp_dist;
				if(hash[tyi*m_width+txi] == 0)
				{
					pq_pixel.push(Pixel(txi, tyi));
					hash[tyi*m_width+txi] = 2;
				}
			}
		}
	}
}

void MSLIC::Matric_DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(
    const unsigned int*                             ubuff,
	const int					width,
	const int					height,
	int*&						klabels,
	const int&					K,//required number of superpixels
    const double&               compactness
	)//weight given to spatial distance
{
    Matric_DoSuperpixelSegmentation_ForGivenSuperpixelSize(ubuff,width,height,klabels,K,compactness);
}

void MSLIC::Matric_DoSuperpixelSegmentation_ForGivenSuperpixelSize(
    const unsigned int*         ubuff,
	const int					width,
	const int					height,
	int*&						klabels,
    const int&					K,
    const double&               compactness)
{
    //------------------------------------------------
	double avg_area = double(width*height)/double(K);
    const int STEP = sqrt(avg_area)+0.5;
	//invwt = 1.0/((STEP/compactness)*(STEP/compactness));
	invwt = compactness*compactness/avg_area;
	//cout<<"invwt: "<<invwt<<endl;
    //------------------------------------------------
	m_width  = width;
	m_height = height;
	sz = m_width*m_height;
	//klabels.resize( sz, -1 );
	//--------------------------------------------------
	klabels = new int[sz];
	for( int s = 0; s < sz; s++ ) klabels[s] = -1;
    //--------------------------------------------------
    if(1)//LAB, the default option
    {
        DoRGBtoLABConversion(ubuff, m_lvec, m_avec, m_bvec);
    }
    else//RGB
    {
        m_lvec = new double[sz]; m_avec = new double[sz]; m_bvec = new double[sz];
        for( int i = 0; i < sz; i++ )
        {
                m_lvec[i] = ubuff[i] >> 16 & 0xff;
                m_avec[i] = ubuff[i] >>  8 & 0xff;
                m_bvec[i] = ubuff[i]       & 0xff;
        }
    }
	//--------------------------------------------------
	Matric_MakeMatric();
	//--------------------------------------------------
    bool perturbseeds(false);
	vector<double> edgemag(0);
	if(perturbseeds) DetectLabEdges(m_lvec, m_avec, m_bvec, m_width, m_height, edgemag);

	//GetLABXYSeeds_ForGivenStepSize(kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, STEP, perturbseeds, edgemag);
	//Matric_GetLABXYSeeds_ForGivenStepSize(STEP);
	Matric_GetLABXYSeeds_ForGivenStepSize_kmeanspp(K);

	int numk=kseedsx.size();
	vector<double> adaptk(numk, 1);
	vector<double> newadaptk(numk, 1);
	for (int itr=0;itr<iteration;itr++)
	{	
		//compute dist & vt
		Matric_PerformSuperpixelMSLIC(klabels, STEP,adaptk,compactness);


		vector<int> klabels_debug;
		for( int s = 0; s < sz; s++ ) klabels_debug.push_back(klabels[s]);
		int numlabels = kseedsx.size();
		//cout<<numlabels<<" ";
		//merge small cell
		//cout<<numlabels<<" ";

		//Matric_EnforceLabelConnectivity(itr, klabels, K);

		//compute center and split
		if(itr != iteration-1)
			Matric_SuperpixelSplit(itr, klabels, adaptk, numlabels);
		else
		{
			Matric_SuperpixelSplit(itr, klabels, adaptk, numlabels);
			Matric_EnforceLabelConnectivity(itr, klabels, K);
		}

		numlabels = kseedsx.size();
		//cout<<numlabels<<endl;
	}
}

void MSLIC::Matric_GetLABXYSeeds_ForGivenStepSize(const int& STEP)
{
	int numseeds(0);
	int n(0);

	int xstrips = (0.5+double(m_width)/double(STEP));
	int ystrips = (0.5+double(m_height)/double(STEP));

    int xerr = m_width  - STEP*xstrips;if(xerr < 0){xstrips--;xerr = m_width - STEP*xstrips;}
    int yerr = m_height - STEP*ystrips;if(yerr < 0){ystrips--;yerr = m_height- STEP*ystrips;}

	double xerrperstrip = double(xerr)/double(xstrips);
	double yerrperstrip = double(yerr)/double(ystrips);

	int xoff = STEP/2;
	int yoff = STEP/2;
	//-------------------------
	numseeds = xstrips*ystrips;
	//-------------------------
	kseedsx.resize(numseeds);
	kseedsy.resize(numseeds);

	for( int y = 0; y < ystrips; y++ )
	{
		int ye = y*yerrperstrip;
		for( int x = 0; x < xstrips; x++ )
		{
			int xe = x*xerrperstrip;
            int seedx = (x*STEP+xoff+xe);
            int seedy = (y*STEP+yoff+ye);
            int i = seedy*m_width + seedx;
			
            kseedsx[n] = seedx;
            kseedsy[n] = seedy;
			n++;
		}
	}
}

void MSLIC::Matric_GetLABXYSeeds_ForGivenStepSize_kmeanspp(const int& K)
{
	double dur;
	clock_t start,end;
	start = clock();
	//code begin
	kseedsx.resize(K);
	kseedsy.resize(K);

	vector<double> tmp_distvec(sz, 1e9);
	vector<int> tmp_Dijkstra_hash(sz, 0);

	distvec = vector<double>(sz, 1e9);

	for(int k=0;k<K;k++)
	{
		double sum = 0;
		for(int i=0;i<sz;i++)
			sum += distvec[i]*distvec[i];
		//cout<<sum<<endl;
		double rand_val = rand()/(RAND_MAX+0.1)*sum;
		//cout<<rand_val<<endl;
		sum=0;
		int rand_rank;
		for(int i=0;i<sz;i++)
		{
			sum += distvec[i]*distvec[i];
			if(sum > rand_val)
			{
				rand_rank = i;
				break;
			}
		}
		//cout<<rand_rank<<endl;
		kseedsx[k] = rand_rank%m_width;
		kseedsy[k] = rand_rank/m_width;
		//cout<<kseedsx[k]<<" "<<kseedsy[k]<<endl;
		int offset;
		offset =50;
		int offsetx = offset;
		int offsety = offset;
		int y1 = max(0.0,			kseedsy[k]-offsety);
		int y2 = min((double)m_height,	kseedsy[k]+offsety);
		int x1 = max(0.0,			kseedsx[k]-offsetx);
		int x2 = min((double)m_width,	kseedsx[k]+offsetx);
		Matric_Dynamic_Dist(tmp_distvec, tmp_Dijkstra_hash, kseedsx[k], kseedsy[k], x1, y1, x2, y2);
		for( int y = y1; y < y2; y++ )
			for( int x = x1; x < x2; x++ )
			{
				int i = y*m_width + x;
				if( tmp_distvec[i] < distvec[i] )
				{
					distvec[i] = tmp_distvec[i];
				}
			}
	}
	//code end
	end = clock();
	dur = (double)(end - start);
	//printf("Use Time:%f\n",(dur/CLOCKS_PER_SEC));
	/*ofstream fout(string("center")+to_string(iteration));
	for(int i=0;i<K;i++)
		fout<<kseedsx[i]<<" "<<kseedsy[i]<<endl;
	fout.close();*/
}

void MSLIC::Matric_PerformSuperpixelMSLIC(
	int*&						klabels,
	const int					STEP,
	vector<double>&				adaptk,
	const double&				M)
{
	const int numk = kseedsx.size();
	//----------------
	int offset = STEP;
	//if(STEP < 8) offset = STEP*1.5;//to prevent a crash due to a very small step size
	//----------------

	//vector<double> distvec(sz, DBL_MAX);
	distvec = vector<double>(sz, DBL_MAX);
	vector<double> tmp_distvec(sz, DBL_MAX);
	vector<int> tmp_Dijkstra_hash(sz, 0);


	for( int n = 0; n < numk; n++ )
	{
		int speed = 2;
		offset=(int)speed*(STEP+5)*adaptk[n];//offset is undetermined

		int y0 = kseedsy[n];
		int x0 = kseedsx[n];
		int y1 = max(0.0,			kseedsy[n]-offset);
		int y2 = min((double)m_height,	kseedsy[n]+offset);
		int x1 = max(0.0,			kseedsx[n]-offset);
		int x2 = min((double)m_width,	kseedsx[n]+offset);

		Matric_Dynamic_Dist_Without_Segment(klabels, tmp_Dijkstra_hash, x0, y0, x1, y1, x2, y2, n);

		/*Matric_Dynamic_Dist(tmp_distvec, tmp_Dijkstra_hash, x0, y0, x1, y1, x2, y2);

		for( int y = y1; y < y2; y++ )
		{
			for( int x = x1; x < x2; x++ )
			{
				int i = y*m_width + x;
				//cout<<tmp_distvec[i]<<endl;
				if( tmp_distvec[i] < distvec[i] )
				{
					distvec[i] = tmp_distvec[i];
					klabels[i]  = n;
				}
			}
		}*/
	}
}

void MSLIC::Matric_EnforceLabelConnectivity(
	int							itr,
	int*&						klabels,
	const int&					K
	) 
{
	vector<int> nlabels(sz, -1);
	vector<int> hash(sz, 0);
	for(int k=0;k<K;k++)
	{
		queue<Pixel> pq_pixel;
		pq_pixel.push(Pixel(kseedsx[k], kseedsy[k]));
		while(!pq_pixel.empty())
		{
			Pixel tmp_pixel = pq_pixel.front();
			int tx = tmp_pixel.x;
			int ty = tmp_pixel.y;
			pq_pixel.pop();
			nlabels[ty*m_width+tx] = k;
			hash[ty*m_width+tx] = 1;
			//go through all neighbor
			for(int i=0;i<8;i++)
			{
				int txi = tx + nb_dx[i];
				int tyi = ty + nb_dy[i];
				if(txi < 0 || txi >= m_width || tyi < 0 || tyi >= m_height)
					continue;
				if(klabels[tyi*m_width+txi] != k)
					continue;
				if(hash[tyi*m_width+txi] == 1)
					continue;
				if(hash[tyi*m_width+txi] == 0)
				{
					pq_pixel.push(Pixel(txi, tyi));
					hash[tyi*m_width+txi] = 2;
				}
			}
		}
	}
	int co = 0;
	for(int i=0;i<sz;i++)
	{
		if(nlabels[i] == -1)
		{
			co++;
			int tmplabel;
			if((i%m_width)>0)
				tmplabel = nlabels[i-1];
			else if((i/m_width)>0)
				tmplabel = nlabels[i-m_width];
			else if(klabels[i] != -1)
				tmplabel = klabels[i];
			else
				tmplabel = 0;
			queue<Pixel> pq_pixel;
			pq_pixel.push(Pixel(i%m_width, i/m_width));
			while(!pq_pixel.empty())
			{
				Pixel tmp_pixel = pq_pixel.front();
				int tx = tmp_pixel.x;
				int ty = tmp_pixel.y;
				//cout<<tx<<" "<<ty<<endl;
				pq_pixel.pop();
				nlabels[ty*m_width+tx] = tmplabel;
				hash[ty*m_width+tx] = 1;
				//go through all neighbor
				for(int i=0;i<8;i++)
				{
					int txi = tx + nb_dx[i];
					int tyi = ty + nb_dy[i];
					if(txi < 0 || txi >= m_width || tyi < 0 || tyi >= m_height)
						continue;
					if(nlabels[tyi*m_width+txi] != -1)
						continue;
					if(hash[tyi*m_width+txi] == 1)
						continue;
					if(hash[tyi*m_width+txi] == 0)
					{
						pq_pixel.push(Pixel(txi, tyi));
						hash[tyi*m_width+txi] = 2;
					}
				}
			}
		}
	}
	if(co>0)
		cout<<"fragment: "<<co<<endl;

}



void MSLIC::Matric_SuperpixelSplit(
	int itr,
	int*&						klabels,
	vector<double>&				adaptk,
	int							numk)
{
	vector<double> clustersizex(numk, 0);
	vector<double> clustersizey(numk, 0);
	vector<double> invx(numk, 0);//to store 1/clustersize[k] values
	vector<double> invy(numk, 0);//to store 1/clustersize[k] values
	vector<double> sigmax(numk, 0);
	vector<double> sigmay(numk, 0);

	{int ind(0);
	for( int r = 0; r < m_height; r++ )
	{
		for( int c = 0; c < m_width; c++ )
		{
			if (klabels[ind]>=0)
			{
				int p0 = r*m_width + c;
				double powerx = 1;//1+distvec[p0];
				double powery = 1;//1+distvec[p0];
				//double powerx = 1/(1+distvec[p0]);
				//double powery = 1/(1+distvec[p0]);
				sigmax[klabels[ind]] += c * powerx;
				sigmay[klabels[ind]] += r * powery;
				clustersizex[klabels[ind]] += powerx;
				clustersizey[klabels[ind]] += powery;
				ind++;
			}
		}
	}}



	{for( int k = 0; k < numk; k++ )
	{
		if( clustersizex[k] <= 0 ) clustersizex[k] = 1;
		invx[k] = 1.0/clustersizex[k];//computing inverse now to multiply, than divide later
		if( clustersizey[k] <= 0 ) clustersizey[k] = 1;
		invy[k] = 1.0/clustersizey[k];//computing inverse now to multiply, than divide later
	}}



	kseedsx.clear();
	kseedsy.clear();
	for (int i=0;i<numk;i++)
	{
		kseedsx.push_back(0);
		kseedsy.push_back(0);
	}
	for( int k = 0; k < numk; k++ )
	{
		kseedsx[k] = int(sigmax[k]*invx[k] + 0.5);
		kseedsy[k] = int(sigmay[k]*invy[k] + 0.5);
	}
}