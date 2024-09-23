gfortran -std=legacy test7.f -L/usr/local/lib -lfmm3d -lfmm3dbie
./a.out
python plot_res_hels_plas.py
