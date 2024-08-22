import numpy as np
import numcodecs
import zarr
import os

input_file='/glade/derecho/scratch/haiyingx/Run_Corona_1728x1024x1024_ASD/3D/result_prim_0.060000'
out_file='/glade/derecho/scratch/haiyingx/result_prim_0.060000.bin'
forder_np_arr=np.fromfile(input_file,dtype=np.float32)
print(forder_np_arr.size)
forder_np_arr_3d=forder_np_arr.reshape([1024,1024,1728]).transpose([2,1,0])
forder_np_arr_3d.tofile(out_file)

os.system("./install/bin/sperr3d -c --ftype 32 --dims 1024 1024 1728 --chunks 256 256 288 --bitstream out.bin --decomp_f out.decomp --psnr 102 /glade/derecho/scratch/haiyingx/Run_Corona_1728x1024x1024_ASD/3D/result_prim_0.060000")
