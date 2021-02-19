# Fast Computation of Content-Sensitive Superpixels and Supervoxels Using Q-Distances
This repository contains the source code for the paper [Fast Computation of Content-Sensitive Superpixels and Supervoxels Using Q-Distances](https://openaccess.thecvf.com/content_ICCV_2019/html/Ye_Fast_Computation_of_Content-Sensitive_Superpixels_and_Supervoxels_Using_Q-Distances_ICCV_2019_paper.html) by [Yong-Jin Liu](https://cg.cs.tsinghua.edu.cn/people/~Yongjin/Yongjin.htm).


## Running Environment
The developing environment we used is Microsoft Visual C++ 2012, and we have tested our executable files on 64-bit MS Win7, Win 8 and Win 10. We recommend users to compile the process in x64 and Release version.


## subprogram
We define many functions for different tasks in ```MainProcess.h```, and the corresponding input parameters are defined in ```main.cpp```. Users can comment and uncomment the code to compile different versions for different tasks.

### Superpixels for RGB images
You can see the function ```superPixel_rgb_iter_mt```, The command is as follows：
```
qd-CSS_rgb.exe  input_folder  output_folder  number_of_superpixels  iter_max
```
Example: 
```
qd-CSS_rgb.exe examples_input examples_output 300 10
```

### Superpixels for RGBD images

### Supervoxels for RGB images

## Additional notes

Citation:

```
@InProceedings{Ye_2019_ICCV,
author = {Ye, Zipeng and Yi, Ran and Yu, Minjing and Liu, Yong-Jin and He, Ying},
title = {Fast Computation of Content-Sensitive Superpixels and Supervoxels Using Q-Distances},
booktitle = {Proceedings of the IEEE/CVF International Conference on Computer Vision (ICCV)},
month = {October},
year = {2019}
}
```