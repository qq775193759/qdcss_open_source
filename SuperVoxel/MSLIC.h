#if !defined(_MSLIC_H_INCLUDED_)
#define _MSLIC_H_INCLUDED_


#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
using namespace std;


class MSLIC  
{
public:
	MSLIC();
	virtual ~MSLIC();
	//============================================================================
	// Superpixel segmentation for a given step size (superpixel size ~= step*step)
	//============================================================================
        void DoSuperpixelSegmentation_ForGivenSuperpixelSize(
        const unsigned int*                            ubuff,//Each 32 bit unsigned int contains ARGB pixel values.
		const int					width,
		const int					height,
		int*&						klabels,
		int&						numlabels,
                const int&					superpixelsize,
                const double&                                   compactness,
				const double&  merge,
				const double&  split,
				const int& speed);
	//============================================================================
	// Superpixel segmentation for a given number of superpixels
	//============================================================================
        void DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(
        const unsigned int*                             ubuff,
		const int					width,
		const int					height,
		int*&						klabels,
		int&						numlabels,
                const int&					K,//required number of superpixels
                const double&                                   compactness,
				const double&  merge,
				const double&  split,
				const int& speed);//10-20 is a good value for CIELAB space
	//============================================================================
	// Save superpixel labels in a text file in raster scan order
	//============================================================================
	void SaveSuperpixelLabels(
		const int*&					labels,
		const int&					width,
		const int&					height,
		const string&				filename,
		const string&				path);
	//============================================================================
	// Function to draw boundaries around superpixels of a given 'color'.
	//============================================================================
	void DrawContoursAroundSegments(
		unsigned int*&				segmentedImage,
		int*&						labels,
		const int&					width,
		const int&					height,
		const unsigned int&			color );
public:
	int  iteration;
	bool ifdraw;
public:
	vector<double> kseedsl;
	vector<double> kseedsa;
	vector<double> kseedsb;
	vector<double> kseedsx;
	vector<double> kseedsy;
private:
	//============================================================================
	// The main MSLIC algorithm for generating superpixels
	//============================================================================
	void printlabel(int*& klabels,int width,int height,int itr);
	void GetKValues_LABXYZ(
		vector<double>&				kseedsl,
		vector<double>&				kseedsa,
		vector<double>&				kseedsb,
		vector<double>&				kseedsx,
		vector<double>&				kseedsy,
		vector<double>&				kseedsz,
		const int&					STEP);
	//============================================================================
	// Move the superpixel seeds to low gradient positions to avoid putting seeds
	// at region boundaries.
	//============================================================================
	void PerturbSeeds(
		vector<double>&				kseedsl,
		vector<double>&				kseedsa,
		vector<double>&				kseedsb,
		vector<double>&				kseedsx,
		vector<double>&				kseedsy,
		const vector<double>&		edges);
	//============================================================================
	// Detect color edges, to help PerturbSeeds()
	//============================================================================
	void DetectLabEdges(
		const double*				lvec,
		const double*				avec,
		const double*				bvec,
		const int&					width,
		const int&					height,
		vector<double>&				edges);
	//============================================================================
	// sRGB to XYZ conversion; helper for RGB2LAB()
	//============================================================================
	void RGB2XYZ(
		const int&					sR,
		const int&					sG,
		const int&					sB,
		double&						X,
		double&						Y,
		double&						Z);
	//============================================================================
	// sRGB to CIELAB conversion (uses RGB2XYZ function)
	//============================================================================
	void RGB2LAB(
		const int&					sR,
		const int&					sG,
		const int&					sB,
		double&						lval,
		double&						aval,
		double&						bval);
	//============================================================================
	// sRGB to CIELAB conversion for 2-D images
	//============================================================================
	void DoRGBtoLABConversion(
		const unsigned int*&		ubuff,
		double*&					lvec,
		double*&					avec,
		double*&					bvec);
	//============================================================================
	// sRGB to CIELAB conversion for 3-D volumes
	//============================================================================
	void DoRGBtoLABConversion(
		unsigned int**&				ubuff,
		double**&					lvec,
		double**&					avec,
		double**&					bvec);
	//============================================================================
	// Post-processing of MSLIC segmentation, to avoid stray labels.
	//============================================================================
	

public:
	//======================================================================================================
	// Matric Space method
	// coded by YZP
	// d(x, y) : straight stretch in parameter space with a anisotropy matric
	//======================================================================================================
		//------------------------------
		//small function
		//------------------------------
		double Matric_Dist(int x1, int y1, int x2, int y2);//for conputing dist
		void Matric_MakeMatric();
		void Matric_Dijkstra_Dist(vector<double> &res, vector<int> &hash, int x0, int y0, int x1, int y1, int x2, int y2);//for conputing dist
		void Matric_Dynamic_Dist(vector<double> &res, vector<int> &hash, int x0, int y0, int x1, int y1, int x2, int y2);//for conputing dist
		void Matric_Dynamic_Dist_Without_Segment(int*& klabels, vector<int> &hash, int x0, int y0, int x1, int y1, int x2, int y2, int label_no);//for conputing dist
	//============================================================================
	// Superpixel segmentation for a given step size (superpixel size ~= step*step)
	//============================================================================
        void Matric_DoSuperpixelSegmentation_ForGivenSuperpixelSize(
        const unsigned int*         ubuff,//Each 32 bit unsigned int contains ARGB pixel values.
		const int					width,
		const int					height,
		int*&						klabels,
        const int&					K,
        const double&               compactness
		);
	//============================================================================
	// Superpixel segmentation for a given number of superpixels
	//============================================================================
        void Matric_DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(
        const unsigned int*                             ubuff,
		const int					width,
		const int					height,
		int*&						klabels,
        const int&					K,//required number of superpixels
        const double&                                   compactness
		);//10-20 is a good value for CIELAB space
	//============================================================================
	// seed spread
	//============================================================================
	void Matric_GetLABXYSeeds_ForGivenStepSize(const int& STEP);
	void Matric_GetLABXYSeeds_ForGivenStepSize_kmeanspp(const int& K);
	//============================================================================
	// Compute VT using matric
	//============================================================================
		void Matric_PerformSuperpixelMSLIC(
		int*&						klabels,
		const int					STEP,
		vector<double>&	    adaptk,
		const double&				m = 10.0
		);
	//============================================================================
	// Post-processing of MSLIC segmentation, to avoid stray labels.
	//============================================================================
		void Matric_EnforceLabelConnectivity(
		int							itr,
		int*&						klabels,
		const int&					K
		); 
	//============================================================================
	// Compute center using matric
	//============================================================================
		void Matric_SuperpixelSplit(
		int itr,
		int*&						klabels,
		vector<double>&				adaptk,
		int							numk
		);

private:
	int										m_width;//x
	int										m_height;//y
	int										m_depth;

	double*									m_lvec;
	double*									m_avec;
	double*									m_bvec;

	double**								m_lvecvec;
	double**								m_avecvec;
	double**								m_bvecvec;

	double invwt;
	int sz;

	//dist
	vector<double> distNeighbor[8];
	int nb_dx[8];
	int nb_dy[8];
	//  7 0 4 
	//  3   1
	//  6 2 5
	//path length
	vector<double> distvec;

public:
	//depth
	int depth_flag;
	vector<int> depth_vec;
	double depth_power;
};

#endif // !defined(_MSLIC_H_INCLUDED_)
