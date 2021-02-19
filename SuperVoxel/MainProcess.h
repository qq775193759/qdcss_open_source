#pragma once

#include <opencv2\opencv.hpp>
#include "MSLIC.h"
#include "SuperVoxel_MatricIMSLIC.h"
#include <iostream>
#include <fstream>
#include <string>


void superPixel_rgb_once(const string pic_name, const string out_name, int m_compactness, int K);
void superPixel_rgb_iter_once(const string pic_name, const string out_name, int K, int iter);
void superPixel_rgbd_autodp_once(const string pic_name, const string out_name, double depth_power, int m_compactness, int K);
void superPixel_rgbd_iter_once(const string pic_name, const string depth_name, const string out_name, int K, int iter);




void superPixel_rgb(const string pic_path, const string out_path, int m_compactness = 20.0, int K = 300);
void superPixel_rgb_iter(const string pic_path, const string out_path, int K, int iter);
void superPixel_rgb_iter_mt(const string pic_path, const string out_path, int K, int iter);
void superPixel_rgbd_autodp(const string pic_path, const string out_path, double depth_power = 1e-3, int m_compactness = 20.0, int K = 300);
void superPixel_rgbd_iter(const string pic_path, const string depth_path, const string out_path, int K, int iter);
void superPixel_rgbd_iter_mt(const string pic_path, const string depth_path, const string out_path, int K, int iter);




void superVoxel_readPath(const string pic_path, SuperVoxel_MatricIMSLIC& imslic);
void superVoxel_rgb(const string pic_path, const string out_path, const int pK, const double pC_xy, const double pC_z, const int iters);
void superVoxel_submit(const string pic_path, const string out_path, const int pK, const int iters);