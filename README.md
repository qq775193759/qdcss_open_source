# Fast Computation of Content-Sensitive Superpixels and Supervoxels Using Q-Distances
This repository contains the source code for the ICCV 2019 paper [Fast Computation of Content-Sensitive Superpixels and Supervoxels Using Q-Distances](https://openaccess.thecvf.com/content_ICCV_2019/html/Ye_Fast_Computation_of_Content-Sensitive_Superpixels_and_Supervoxels_Using_Q-Distances_ICCV_2019_paper.html) by [Yong-Jin Liu](https://cg.cs.tsinghua.edu.cn/people/~Yongjin/Yongjin.htm).


## Running Environment
The developing environment we used is Microsoft Visual C++ 2012, and we have tested our executable files on 64-bit MS Win7, Win 8 and Win 10. We recommend users to compile the process in x64 and Release version.


## subprogram
We define many functions for different tasks in ```MainProcess.h```, and the corresponding input parameters are defined in ```main.cpp```. Users can comment and uncomment the code to compile different versions for different tasks.

### Superpixels for RGB images
You can see the function ```superPixel_rgb_iter_mt```, The command is as follows (rename the executable file as ```qd-CSS_rgb.exe```):
```
qd-CSS_rgb.exe  input_folder  output_folder  number_of_superpixels  iter_max
```
input_folder: the filepath of input folder containing RGB images which would be segmented.The types can be JPG, PNG or BMP.
output_folder: the folder filepath contains corresponding superpixels result BMP images.
number_of_superpixels: the number of superpixels you want to generate, which should be integer. 200-700 recommended.
iter_max: the max number of interation, which should be integer. 5-20 recommended.

Example: 
```
qd-CSS_rgb.exe examples_input examples_output 300 10
```

### Superpixels for RGBD images
You can see the function ```superPixel_rgbd_iter_mt```, The command is as follows (rename the executable file as ```qd-CSS_rgbd.exe```):
```
qd-CSS_rgbd.exe  rgb_image_folder  depth_image_folder  output_folder  number_of_superpixels  iter_max
```
rgb_image_fold: the filepath of input folder containing RGB images which would be segmented.The types should be PNG format.
depth_image_folder: the folder filepath contains corresponding depth input which should also be PNG format.
output_folder: the folder filepath contains corresponding superpixels result BMP images.
number_of_superpixels: the number of superpixels you want to generate, which should be integer. 200-700 recommended.
iter_max: the max number of interation, which should be integer. 5-20 recommended.

Example: 
```
qd-CSS_rgbd.exe images depth output 300 10
```

### Supervoxels for RGB videos
You can see the function ```superVoxel_submit```, The command is as follows (rename the executable file as ```qd-CSS_sv.exe```):
```
qd-CSS_sv.exe  input_folder  output_folder  number_of_supervoxels  iter_max
```
input_folder: the filepath of input folder containing video frames which would be segmented.The types can be JPG, PNG or BMP.
output_folder: the folder filepath contains corresponding supervoxels result png images.
number_of_supervoxels: the number of supervoxels you want to generate, which should be integer. 200-5000 recommended.
iter_max: the max number of interation, which should be integer. 5-20 recommended.

Example: 
```
qd-CSS_sv.exe sv_input sv_output 500 10
```

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