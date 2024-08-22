import os
import numpy as np
import time
import zarr

#start this script:
#module load conda
#conda activate derecho_sperr103


enable_sperr=0

var_name_list = [
    "yz_slice_0132",
    "xz_slice_0512",
    "xy_slice_0512",
    "yz_slice_0148",
    "yz_slice_0180",
    "yz_slice_0244",
    #"tau_slice_0.000001",
    #"tau_slice_0.000100",
    #"tau_slice_0.001",
    #"tau_slice_0.010",
    #"tau_slice_0.100",
    #"tau_slice_1.000",
    "yz_slice_0000",
    "yz_slice_0058",
    "yz_slice_0116",
    "yz_slice_0372",
    "yz_slice_0628",
    "yz_slice_1140",
    "yz_slice_1652",
    #"I_out"
]
iter_list=['060100','060200','060300','060400']
iter_list=['060400']
dir_name='/glade/derecho/scratch/haiyingx/Run_Corona_1728x1024x1024_ASD'
orig_in_path=f'{dir_name}/2D/'
zarr_in_path=f'{dir_name}/2D/'
scratch=os.environ.get('SCRATCH')
sperr_out_path=f'{scratch}/sperr_data/'
zarr_filename=f'{zarr_in_path}/result_prim.4.lossless.single.122.zarr'
#zarr_filename=f'{zarr_in_path}{var_name}.n4.zarr'

def sperr_compress(filename,in_path,out_path):
    input_file=f'{in_path}/{filename}'
    input_file_corder=f'{input_file}.c_order'
    comp_file=f'{out_path}/{filename}.c_order.bin'
    out_file=f'{out_path}/{filename}.c_order.decomp'
    forder_np_arr=np.fromfile(input_file,dtype=np.float32)
    print(forder_np_arr.size)
    forder_np_arr_3d=forder_np_arr.reshape([1024,1024,1728]).transpose([2,1,0])
    forder_np_arr_3d.tofile(f'{input_file}.c_order')
    print(input_file_corder)

    #cmd='export LD_LIBRARY_PATH=/glade/work/haiyingx/conda-envs/derecho_sperr103/lib/python3.9/site-packages/numcodecs/../lib64:$LD_LIBRARY_PATH'
    os.system(cmd)
    cmd='/glade/derecho/scratch/haiyingx/numcodecs_sperr_derecho/SPERR/install/bin/sperr3d -c --ftype 32 --dims 1024 1024 1728 --chunks 256 256 288 --bitstream {} --decomp_f {} --psnr 130 {}'.format(comp_file,out_file,input_file_corder)
    print(cmd)
    start = time.perf_counter()
    os.system(cmd)
    end = time.perf_counter()
    print("Sperr compression time = ",end-start)
    #os.system('rm {}'.format(comp_file))

def zarr_sperr_compress(var_name,niter,in_path,out_path):
    input_file=f'{in_path}/{var_name}.{niter}'

    shape1=(1728,1024,1024)
    chunks1=(288,256,256)

    dim3D=['1728','1024','1024']
    print('input_file=',input_file)
    f1=np.fromfile(input_file,dtype=np.float32)

    f1_3d=f1.reshape([1024,1024,1728]).transpose([2,1,0])
    print(f1_3d[0,0,0:10])
    zarr_compression_kw = dict(compressor=zarr.Sperr(mode=2,level=130,pre=''))
    zarr_output_filename=f'{out_path}/{var_name}.{niter}.zarr'
    start = time.perf_counter()
    store = zarr.DirectoryStore(zarr_output_filename)

    print('zarr_output_filename=',zarr_output_filename)
    prim=zarr.group(store=store,overwrite=True)
    if enable_sperr:
        prim_0=prim.create_dataset(var_name,shape=shape1,chunks=chunks1,dtype='f4',**zarr_compression_kw,order='C')
    else:
        prim_0=prim.create_dataset(var_name,shape=shape1,chunks=chunks1,dtype='f4')
    prim_0.attrs['_ARRAY_DIMENSIONS']=dim3D
    prim_0[:,:,:]=f1_3d
    end = time.perf_counter()
    print("Done, zarr compression time = ",end-start)


    #f1_read=np.fromfile(fpath+'/3D/'+'sperr.result_prim_0.bin',dtype=np.float32)
    #f1_read_3d=f1_read.reshape([1024,1024,1728]).transpose([2,1,0])
    #zarr_read=zarr.open(zarr_filename)
    #zarr_read_3d=zarr_read[var_name]
    #print(z1_read[10,0,0:10])
    #print(np.array_equal(zarr_read_3d,f1_3d))

def compare_orig_zarr(orig_filename,zarr_filename,var_name):
    print(orig_filename)
    print(zarr_filename,var_name)
    orig_read=np.fromfile(orig_filename,dtype=np.float32)
    print(orig_read[0:4])
    if var_name.startswith('xy'): 
        orig_read_3d=orig_read[4:].reshape([10,1728,1024])#.transpose([2,1,0])
    elif var_name.startswith('xz'):
        orig_read_3d=orig_read[4:].reshape([10,1728,1024])
    elif var_name.startswith('yz'):
        orig_read_3d=orig_read[4:].reshape([1024,1024,10]).transpose([2,1,0])
    elif var_name.startswith('tau'):
        orig_read_3d=orig_read[4:].reshape([1024,1024])#.transpose([2,1,0])
    print('orig ',orig_read_3d[0,256,1014:1024])
    start = time.perf_counter()
    zarr_read=zarr.open(zarr_filename)
    zarr_read_3d=zarr_read[var_name]
    print('zarr ',zarr_read_3d[0,256,1014:1024])
    end = time.perf_counter()
    print("Done, zarr decompression time = ",end-start)
    print(np.array_equal(zarr_read_3d,orig_read_3d))

def compare_sperr_zarr(sperr_filename,zarr_filename,var_name):
    print('sperr_filename=',sperr_filename)
    sperr_input_file=np.fromfile(sperr_filename,dtype=np.float32)
    #Fortran array need to be tranposed
    #sperr_arr_3d=sperr_input_file.reshape([1024,1024,1728]).transpose([2,1,0])
    sperr_arr_3d=sperr_input_file.reshape([1728,1024,1024])
    print(sperr_arr_3d[0,0,0:10])
    print('zarr_filename=',zarr_filename)
    zarr_input_file=zarr.open(zarr_filename)
    zarr_arr_3d=zarr_input_file[var_name]
    print(zarr_arr_3d[0,0,0:10])
    print(np.array_equal(sperr_arr_3d, zarr_arr_3d))
   
def test_adios2io_correctness():
  ni=3
  for i in iter_list:
     ni=ni+1
     for n in var_name_list:
        if not enable_sperr:
            if n.startswith('tau'):
                compare_orig_zarr(f'{orig_in_path}/{n}.{i}',f'{zarr_in_path}/{n}_0.{ni}.lossless.single.122.zarr', f'{n}')
            else:
                compare_orig_zarr(f'{orig_in_path}/{n}.{i}',f'{zarr_in_path}/{n[:8]}.{ni}.lossless.single.122.zarr', f'{n}')


def test_commandline_sperr():
  for i in iter_list:
    for n in var_name_list:
       src_name=f'{n}.{i}'
       sperr_compress(src_name,orig_in_path,sperr_out_path)

def test_compare_sperr_zarr():
  ni=0
  for i in iter_list:
   if ni==0:
    ni=ni+1
    for n in var_name_list:
      generated_zarr_filename=f'{zarr_in_path}/{n}.{i}.zarr'
      sperr_file=f'{sperr_out_path}/{n}.{i}.c_order.decomp'
      compare_sperr_zarr(sperr_file,generated_zarr_filename,n)

def test_zarr_sperr_compression():
  ni=0
  for i in iter_list:
    ni=ni+1
    if ni==1:
      for n in var_name_list:
        print(n)
        zarr_sperr_compress(n,i,orig_in_path,zarr_in_path)
   
test_adios2io_correctness()
#test_zarr_sperr_compression()
#test_commandline_sperr()
#test_compare_sperr_zarr()
