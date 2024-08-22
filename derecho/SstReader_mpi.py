from dask_mpi import initialize
import dask.array as da
from distributed import Client, wait
import json
from operator import floordiv
from mpi4py import MPI
import numpy as np
import adios2
import zarr
from parse_parameters import parsed_param_dict,parsed_backup_dict,get_output_varname,format_level


with open('slice_ranks.txt', 'r') as srjson:
    rank_info = json.load(srjson)
print(rank_info)
begin_iter='000000'
#begin_globiter=60000
need_zip=0
begin_globiter=int(parsed_backup_dict(""))
param_dict=parsed_param_dict("")
multifile = 1
singlefile = 0
loop_order = np.array([[ 1, 2, 0 ],[ 0, 1, 2 ],[ 2, 1, 0]])
print(param_dict['maxiter'],param_dict['resfreq'])
nsteps = (param_dict['maxiter'][0]-begin_globiter)//param_dict['resfreq'][0]
print("nsteps ",nsteps)
enable_double = param_dict['enable_adios2_double'][0] 
print('enable',enable_double,type(enable_double))
procs = param_dict['procs']
gsize = param_dict['gsize']
nslvar=sum(param_dict['xy_var'])
nslvar_yz=sum(param_dict['yz_var'])
print('nslvar=',nslvar)
tau_nslvar=sum(param_dict['tau_var'])
eos_output=param_dict['eos_var']
diag_output=get_output_varname(param_dict['diag_var'],param_dict['diag_output'])
tau_lev=format_level(param_dict['tau_lev'])

xy_lev=param_dict['xy_lev']
format_xy_lev=[]
for i in xy_lev:
    format_xy_lev.append( "{:04d}".format(i))
yz_lev=param_dict['yz_lev']
if len(yz_lev) == 1 and yz_lev[0]==0:
    yz_slice=0
else:
    yz_slice=1
format_yz_lev=[]
for i in yz_lev:
    format_yz_lev.append( "{:04d}".format(i))
dem=param_dict['DEM'][0]
print('dem ',dem,type(dem))

#print('nslvar',nslvar)
nslvar_corona=7
lsize=list(map(floordiv, gsize, procs))
print('lsize',lsize)

d1=loop_order[:,0]
d2=loop_order[:,1]
d3=loop_order[:,2]

# MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()
cartesian3d =comm.Create_cart(dims = procs, periods =[False,False,False],reorder=False)
lrank = cartesian3d.Get_coords(rank)
print('nprocs ',nprocs,'rank ',rank,lrank)

ghosts=[0,0,0]
gbeg=[0,0,0]
start=[]
end=[]
count=[]
for i in range(3):
    start.append(lrank[i]*lsize[i])
    end.append(start[i]+lsize[i])
    count.append(lsize[i])
    
#print(rank,' start ',start)
#print(rank, ' end ',end)
#print(rank, ' count ',count)

#split comm to same as MURaM comm
xy_comm=comm.Split(color=lrank[2],key=rank)
xz_comm=comm.Split(color=lrank[1],key=rank)
yz_comm=comm.Split(color=lrank[0],key=rank)

x_color = lrank[1]+lrank[2]*procs[1]
y_color = lrank[0]+lrank[2]*procs[0]
z_color = lrank[0]+lrank[1]*procs[0]
print("Rank", rank,"x_color",x_color,"y_color",y_color,"z_color",z_color)
xcol_comm = comm.Split(color=x_color,key=rank)
ycol_comm = comm.Split(color=y_color,key=rank)
zcol_comm = comm.Split(color=z_color,key=rank)

xy_rank=xy_comm.Get_rank()
xz_rank=xz_comm.Get_rank()
yz_rank=yz_comm.Get_rank()

yz_procs=yz_comm.Get_size()

xcol_rank=xcol_comm.Get_rank()
ycol_rank=ycol_comm.Get_rank()
zcol_rank=zcol_comm.Get_rank()
print("rank",rank, "xy",xy_rank,"xz",xz_rank,"yz",yz_rank,"xcol",xcol_rank,"ycol",ycol_rank,"zcol",zcol_rank)



# Zarr file compressor definition 
zarr_compression_type, zarr_compression_level = 'gzip', 4
zarr_compression_kw = dict(compressor=zarr.Sperr(mode=2,level=130))
zarr_compression_kw = dict(compressor=zarr.GZip(level=4))
dtype ='f4'
if enable_double:
    dtype='f8'
#print('dtype=',dtype)
chunk_size1 = np.dtype(dtype).itemsize * (lsize[0]*lsize[1]*lsize[2])

# Get the group of processes in MPI_COMM_WORLD
world_group=comm.Get_group()


# Construct a group containing all of the worker ranks in world_group
#worker_group=world_group.Excl([0,1])
worker_group=world_group.Dup()

# Create a new communicator based on the group
worker_group=world_group.Dup()

# Create a new communicator based on the group
worker_comm=comm.Create(worker_group)

# ADIOS MPI Communicator, debug mode
#adios = adios2.ADIOS(worker_comm, adios2.DebugON)
adios = adios2.ADIOS(worker_comm)

# ADIOS IO
bpIO = adios.DeclareIO("io_xyzl")
# ADIOS Engine
bpIO.SetEngine('Sst')
#bpIO.SetEngine('BPFile')
bpIOParams = {}
bpIOParams['Threads'] = '1'
#bpIOParams['ProfileUnits'] = 'Microseconds'
bpIOParams['InitialBufferSize'] = '2GB'
bpIOParams['OpenTimeoutSecs'] = '60000'
#bpIOParams["MarshalMethod"]="FFS"
#bpIOParams["SstVerbose"]= "5"
#bpIOParams["DataTransport"]="RDMA"
bpIOParams["DataTransport"]="WAN"
#bpIOParams["ControlTransport"]="enet"
#bpIOParams["NetworkInterface"]="hsn0"
bpIO.SetParameters(bpIOParams)


bpFileReader = bpIO.Open("result_prim", adios2.Mode.Read)
for i in range(nsteps):
    print("start nsteps ", i)
    ni=i+1
    if multifile:
        dim3D = ['grid0','grid1','grid2']
        dim3D_2 = [ nslvar,'grid1','grid2']
        dim3D_3 = [ nslvar_corona,'grid1','grid2']
        dim2D_1 = [ 'grid1','grid2']
        shape_scalar = (1)
        chunk_scalar = (1)
        shape1 = tuple(gsize) 
        chunks1 = tuple(lsize)
        shape_yz = (nslvar_yz, gsize[1], gsize[2])
        chunks_yz = (nslvar_yz, lsize[1], lsize[2])
        shape3 = (nslvar,gsize[0],gsize[1])
        chunks3 = (nslvar,lsize[0],lsize[1])
        shape4 = (gsize[1],gsize[2])
        chunks4 = (lsize[1],lsize[2])
    else:
        dim3D = ['timestep','grid0','grid1','grid2']
        dim2D = ['timestep','grid1','grid2']
        shape1 = (nsteps, 192, 64, 64)
        chunks1 = (1, 96, 32, 64)
        #shape2 = (nsteps, 8, 64, 64)
        #chunks2 = (1, 8, 32, 64)
    shape2=[]
    chunks2=[]
    for i in range(loop_order[:,0].size):
        d1=loop_order[i][0];
        d2=loop_order[i][1];
        d3=loop_order[i][2];
        chunks2.append([nslvar_corona,lsize[d2],lsize[d3]])
        shape2.append([nslvar_corona,gsize[d2],gsize[d3]])
        print('rank=',rank,'shape2=',shape2[i],chunks2[i])
    if singlefile or multifile:
        store = zarr.DirectoryStore(f'3D/result_prim.{ni}.lossless.single.122.zarr')
        store_I = zarr.DirectoryStore(f'2D/I_out.{ni}.lossless.single.122.zarr')
        store_xy = zarr.DirectoryStore(f'2D/xy_slice.{ni}.lossless.single.122.zarr')
        
        if yz_slice:
            store_yz = zarr.DirectoryStore(f'2D/yz_slice.{ni}.lossless.single.122.zarr')
        store_tau = zarr.DirectoryStore(f'2D/tau_slice.{ni}.lossless.single.122.zarr')
        if dem:
            store_corona = zarr.DirectoryStore(f'2D/corona_emission.{ni}.lossless.single.122.zarr')
        
        if need_zip:
            store = zarr.ZipStore(f'3D/result_prim.{ni}.lossless.single.122.zip', mode='w')
            store_I = zarr.ZipStore(f'2D/I_out.{ni}.lossless.single.122.zip', mode='w')
            store_xy = zarr.ZipStore(f'2D/xy_slice.{ni}.lossless.single.122.zip', mode='w')
            if yz_slice:
                store_yz = zarr.ZipStore(f'2D/yz_slice.{ni}.lossless.single.122.zip', mode='w')
            store_tau = zarr.ZipStore(f'2D/tau_slice.{ni}.lossless.single.122.zip', mode='w')
            if dem:
                store_corona = zarr.ZipStore(f'2D/corona_emission.{ni}.lossless.single.122.zip', mode='w')
                
        if rank == 0:
            prim = zarr.group(store=store, overwrite=True)
            prim_I = zarr.group(store=store_I, overwrite=True)
            prim_xy = zarr.group(store=store_xy, overwrite=True)
            if yz_slice:
                prim_yz = zarr.group(store=store_yz, overwrite=True)
            prim_tau = zarr.group(store=store_tau, overwrite=True)
            if dem:
                prim_corona = zarr.group(store=store_corona, overwrite=True)
            prim_0 = prim.create_dataset('prim_0', shape=shape1, chunks=chunks1,dtype=dtype,**zarr_compression_kw,order='F')
            #prim_0 = prim.create_dataset('prim_0', shape=shape1, chunks=chunks1,dtype=dtype)
            prim_1 = prim.create_dataset('prim_1', shape=shape1, chunks=chunks1,dtype=dtype,**zarr_compression_kw,order='F')
            prim_2 = prim.create_dataset('prim_2', shape=shape1, chunks=chunks1,dtype=dtype,**zarr_compression_kw,order='F')
            prim_3 = prim.create_dataset('prim_3', shape=shape1, chunks=chunks1,dtype=dtype,**zarr_compression_kw,order='F')
            prim_4 = prim.create_dataset('prim_4', shape=shape1, chunks=chunks1,dtype=dtype,**zarr_compression_kw,order='F')
            prim_5 = prim.create_dataset('prim_5', shape=shape1, chunks=chunks1,dtype=dtype,**zarr_compression_kw,order='F')
            prim_6 = prim.create_dataset('prim_6', shape=shape1, chunks=chunks1,dtype=dtype,**zarr_compression_kw,order='F')
            prim_7 = prim.create_dataset('prim_7', shape=shape1, chunks=chunks1,dtype=dtype,**zarr_compression_kw,order='F')
            prim_8 = prim.create_dataset('prim_8', shape=shape1, chunks=chunks1,dtype=dtype,**zarr_compression_kw,order='F')
            diag=[]
            for vi in diag_output:
                diag.append(prim.create_dataset(vi, shape=shape1, chunks=chunks1,dtype=dtype,**zarr_compression_kw,order='F'))
            #diag_0 = prim.create_dataset('Qres', shape=shape1, chunks=chunks1,dtype=dtype,**zarr_compression_kw,order='F')
            #diag_1 = prim.create_dataset('Qvis', shape=shape1, chunks=chunks1,dtype=dtype,**zarr_compression_kw,order='F')
            #diag_2 = prim.create_dataset('Qamb', shape=shape1, chunks=chunks1,dtype=dtype,**zarr_compression_kw,order='F')
            eos=[]
            for vi in eos_output: 
                eos.append(prim.create_dataset(vi, shape=shape1, chunks=chunks1,dtype=dtype,**zarr_compression_kw,order='F'))
            #eos_1 = prim.create_dataset('eosP', shape=shape1, chunks=chunks1,dtype=dtype,**zarr_compression_kw,order='F')
            #eos_2 = prim.create_dataset('eosne', shape=shape1, chunks=chunks1,dtype=dtype,**zarr_compression_kw,order='F')
            #eos_3 = prim.create_dataset('eosrhoi', shape=shape1, chunks=chunks1,dtype=dtype,**zarr_compression_kw,order='F')
            #eos_4 = prim.create_dataset('eosamb', shape=shape1, chunks=chunks1,dtype=dtype,**zarr_compression_kw,order='F')
            #eos_5 = prim.create_dataset('Qtot', shape=shape1, chunks=chunks1,dtype=dtype,**zarr_compression_kw,order='F')
            #eos_6 = prim.create_dataset('tau', shape=shape1, chunks=chunks1,dtype=dtype,**zarr_compression_kw,order='F')
            #eos_7 = prim.create_dataset('Jtot', shape=shape1, chunks=chunks1,dtype=dtype,**zarr_compression_kw,order='F')
            #eos_8 = prim.create_dataset('Stot', shape=shape1, chunks=chunks1,dtype=dtype,**zarr_compression_kw,order='F')
            #eos_9 = prim.create_dataset('QxCor', shape=shape1, chunks=chunks1,dtype=dtype,**zarr_compression_kw,order='F')
           
            for li in tau_lev:
                for vi in range(tau_nslvar):
                    tau_0 = prim_tau.create_dataset(f'tau_slice_{li}_{vi}', shape=shape4, chunks=chunks4,dtype=dtype,**zarr_compression_kw,order='F')
                #tau_1 = prim_tau.create_dataset(f'tau_slice_0100_{vi}', shape=shape4, chunks=chunks4,dtype=dtype,**zarr_compression_kw,order='F')
                #tau_2 = prim_tau.create_dataset(f'tau_slice_0010_{vi}', shape=shape4, chunks=chunks4,dtype=dtype,**zarr_compression_kw,order='F')
            for vi in format_xy_lev:
                xy_slice = prim_xy.create_dataset(f'xy_slice_{vi}', shape=shape3, chunks=chunks3,dtype=dtype,**zarr_compression_kw,order='F')
            if yz_slice:
                for vi in format_yz_lev:
                    yz_slice = prim_yz.create_dataset(f'yz_slice_{vi}', shape=shape_yz, chunks=chunks_yz,dtype=dtype,**zarr_compression_kw,order='F')
            I_out = prim_I.create_dataset('I_out', shape=shape4, chunks=chunks4,dtype=dtype,**zarr_compression_kw,order='F')
            if dem:
                corona_1 = prim_corona.create_dataset('corona_x1',shape=shape2[1],chunks=chunks2[1],dtype=dtype,**zarr_compression_kw)
                corona_2 = prim_corona.create_dataset('corona_x2',shape=shape2[1],chunks=chunks2[1],dtype=dtype,**zarr_compression_kw)
                corona_3 = prim_corona.create_dataset('corona_x3',shape=shape2[1],chunks=chunks2[1],dtype=dtype,**zarr_compression_kw)
                corona_4 = prim_corona.create_dataset('corona_x4',shape=shape2[1],chunks=chunks2[1],dtype=dtype,**zarr_compression_kw)
                corona_y1 = prim_corona.create_dataset('corona_y1',shape=shape2[0],chunks=chunks2[0],dtype=dtype,**zarr_compression_kw)
                corona_y2 = prim_corona.create_dataset('corona_y2',shape=shape2[0],chunks=chunks2[0],dtype=dtype,**zarr_compression_kw)
                corona_y3 = prim_corona.create_dataset('corona_y3',shape=shape2[0],chunks=chunks2[0],dtype=dtype,**zarr_compression_kw)
                corona_y4 = prim_corona.create_dataset('corona_y4',shape=shape2[0],chunks=chunks2[0],dtype=dtype,**zarr_compression_kw)
                corona_z1 = prim_corona.create_dataset('corona_z1',shape=shape2[2],chunks=chunks2[2],dtype=dtype,**zarr_compression_kw)
                corona_z2 = prim_corona.create_dataset('corona_z2',shape=shape2[2],chunks=chunks2[2],dtype=dtype,**zarr_compression_kw)
                corona_z3 = prim_corona.create_dataset('corona_z3',shape=shape2[2],chunks=chunks2[2],dtype=dtype,**zarr_compression_kw)
                corona_z4 = prim_corona.create_dataset('corona_z4',shape=shape2[2],chunks=chunks2[2],dtype=dtype,**zarr_compression_kw)
            prim_0.attrs['_ARRAY_DIMENSIONS']=dim3D
            prim_1.attrs['_ARRAY_DIMENSIONS']=dim3D
            prim_2.attrs['_ARRAY_DIMENSIONS']=dim3D
            prim_3.attrs['_ARRAY_DIMENSIONS']=dim3D
            prim_4.attrs['_ARRAY_DIMENSIONS']=dim3D
            prim_5.attrs['_ARRAY_DIMENSIONS']=dim3D
            prim_6.attrs['_ARRAY_DIMENSIONS']=dim3D
            prim_7.attrs['_ARRAY_DIMENSIONS']=dim3D
            prim_8.attrs['_ARRAY_DIMENSIONS']=dim3D
            for c,vi in enumerate(diag_output):
                diag[c].attrs['_ARRAY_DIMENSIONS']=dim3D
            #diag_1.attrs['_ARRAY_DIMENSIONS']=dim3D
            #diag_2.attrs['_ARRAY_DIMENSIONS']=dim3D
            for ai in range(len(eos_output)):
                eos[ai].attrs['_ARRAY_DIMENSIONS']=dim3D
            for li in tau_lev:
                for vi in range(tau_nslvar):
                    prim_tau[f'tau_slice_{li}_{vi}'].attrs['_ARRAY_DIMENSIONS']=dim2D_1
                #prim_tau[f'tau_slice_0100_{vi}'].attrs['_ARRAY_DIMENSIONS']=dim2D_1
                #prim_tau[f'tau_slice_0010_{vi}'].attrs['_ARRAY_DIMENSIONS']=dim2D_1
            I_out.attrs['_ARRAY_DIMENSIONS']=dim2D_1
            for vi in format_xy_lev:
                prim_xy[f'xy_slice_{vi}'].attrs['_ARRAY_DIMENSIONS']=dim3D_2
            if yz_slice:
                for vi in format_yz_lev:
                    prim_yz[f'yz_slice_{vi}'].attrs['_ARRAY_DIMENSIONS']=dim3D_2
            if dem:
                corona_1.attrs['_ARRAY_DIMENSIONS']=dim3D_3
                corona_2.attrs['_ARRAY_DIMENSIONS']=dim3D_3
                corona_3.attrs['_ARRAY_DIMENSIONS']=dim3D_3
                corona_4.attrs['_ARRAY_DIMENSIONS']=dim3D_3
                corona_y1.attrs['_ARRAY_DIMENSIONS']=dim3D_3
                corona_y2.attrs['_ARRAY_DIMENSIONS']=dim3D_3
                corona_y3.attrs['_ARRAY_DIMENSIONS']=dim3D_3
                corona_y4.attrs['_ARRAY_DIMENSIONS']=dim3D_3
                corona_z1.attrs['_ARRAY_DIMENSIONS']=dim3D_3
                corona_z2.attrs['_ARRAY_DIMENSIONS']=dim3D_3
                corona_z3.attrs['_ARRAY_DIMENSIONS']=dim3D_3
                corona_z4.attrs['_ARRAY_DIMENSIONS']=dim3D_3
        worker_comm.Barrier()
        prim2 = zarr.open_group(store=store,mode='a')
        prim2_I = zarr.open_group(store=store_I,mode='a')
        prim2_xy = zarr.open_group(store=store_xy,mode='a')
        if yz_slice:
            prim2_yz = zarr.open_group(store=store_yz,mode='a')
        prim2_tau = zarr.open_group(store=store_tau,mode='a')
        if dem:
            prim2_corona = zarr.open_group(store=store_corona,mode='a')
        singlefile = 0

    count1 = count
    start1=start
    end1=end
    #print('start ',start)
    #print('count ',count)
    #print('end ',end)

    start2=[]
    count2=[]
    end2=[]
    for i in range(loop_order[:,0].size):
        d1=loop_order[i][0];
        d2=loop_order[i][1];
        d3=loop_order[i][2];
        start2.append([0,lrank[d2]*lsize[d2],
             lrank[d3]*lsize[d3]])
        count2.append([nslvar_corona,lsize[d2],lsize[d3]])
        end2.append([nslvar_corona,lrank[d2]*lsize[d2]+lsize[d2],
             lrank[d3]*lsize[d3]+lsize[d3]])
        print('rank=',rank,'start2=',start2[i],count2[i],end2[i])

    count3=[nslvar,count[0],count[1]] 
    start3=[0,start[0],start[1]]
    end3=[count3[0],end[0],end[1]]
    count_yz=[nslvar_yz,count[1],count[2]] 
    start_yz=[0,start[1],start[2]]
    end_yz=[count_yz[0],end[1],end[2]]
    #print('start3 ',start3,lrank)
    #print('count3 ',count3)
    #print('end3 ',end3)

    count4=[count[1],count[2]]
    #start4=[start2[1],start2[2]]
    start4=[start[1],start[2]]
    end4=[end[1],end[2]]
    #print('start4 ',start4,lrank)
    #print('count4 ',count4)
    #print('end4 ',end4)

    #result_prim data
    if bpFileReader.BeginStep()==adios2.StepStatus.OK:

        var_get = bpIO.InquireVariable(f"prim_0.{begin_iter}")
        var_get1 = bpIO.InquireVariable(f"prim_1.{begin_iter}")
        var_get2 = bpIO.InquireVariable(f"prim_2.{begin_iter}")
        var_get3 = bpIO.InquireVariable(f"prim_3.{begin_iter}")
        var_get4 = bpIO.InquireVariable(f"prim_4.{begin_iter}")
        var_get5 = bpIO.InquireVariable(f"prim_5.{begin_iter}")
        var_get6 = bpIO.InquireVariable(f"prim_6.{begin_iter}")
        var_get7 = bpIO.InquireVariable(f"prim_7.{begin_iter}")
        var_get8 = bpIO.InquireVariable(f"prim_8.{begin_iter}")
        var_get.SetSelection([start,count])
        var_get1.SetSelection([start,count])
        var_get2.SetSelection([start,count])
        var_get3.SetSelection([start,count])
        var_get4.SetSelection([start,count])
        var_get5.SetSelection([start,count])
        var_get6.SetSelection([start,count])
        var_get7.SetSelection([start,count])
        var_get8.SetSelection([start,count])
        inVar=np.zeros(count,dtype=dtype,order='F')
        inVar1=np.zeros(count,dtype=dtype,order='F')
        inVar2=np.zeros(count,dtype=dtype,order='F')
        inVar3=np.zeros(count,dtype=dtype,order='F')
        inVar4=np.zeros(count,dtype=dtype,order='F')
        inVar5=np.zeros(count,dtype=dtype,order='F')
        inVar6=np.zeros(count,dtype=dtype,order='F')
        inVar7=np.zeros(count,dtype=dtype,order='F')
        inVar8=np.zeros(count,dtype=dtype,order='F')
        bpFileReader.Get(var_get, inVar, adios2.Mode.Deferred)
        bpFileReader.Get(var_get1, inVar1, adios2.Mode.Deferred)
        bpFileReader.Get(var_get2, inVar2, adios2.Mode.Deferred)
        bpFileReader.Get(var_get3, inVar3, adios2.Mode.Deferred)
        bpFileReader.Get(var_get4, inVar4, adios2.Mode.Deferred)
        bpFileReader.Get(var_get5, inVar5, adios2.Mode.Deferred)
        bpFileReader.Get(var_get6, inVar6, adios2.Mode.Deferred)
        bpFileReader.Get(var_get7, inVar7, adios2.Mode.Deferred)
        err=bpFileReader.Get(var_get8, inVar8, adios2.Mode.Deferred)
        if err :
            print("var_get8 Get does not works",err)
    
    bpFileReader.EndStep()
    #bpFileReader.PerformGets()

    #eos data
    bpFileReader.BeginStep()
    var_info = bpIO.AvailableVariables()
    for name, info in var_info.items():
        if rank == 2:
            print("attribute_name: " + name)
            for key, value in info.items():
                print("\t" + key + ": " + value)
            print("\n")
    eos_get=[]
    eosVar=[]
    for vv in eos_output:
        print('eos_begin',f'{vv}.{begin_iter}')
        eos_get.append(bpIO.InquireVariable(f'{vv}.{begin_iter}'))
    for vi in range(len(eos_output)):
        eos_get[vi].SetSelection([start,count])
    for vi in range(len(eos_output)):
        eosVar.append(np.zeros(count,dtype=dtype,order='F'))

    for vi in range(len(eos_output)):
        bpFileReader.Get(eos_get[vi], eosVar[vi], adios2.Mode.Sync)
    bpFileReader.EndStep()

    #Diag variables
    bpFileReader.BeginStep()
    diag_get=[]
    diagVar=[]
    for vi in diag_output:
        diag_get.append(bpIO.InquireVariable(f"{vi}.{begin_iter}"))
    for vi,vn in enumerate(diag_output):
        diag_get[vi].SetSelection([start,count])
    for vi,vn in enumerate(diag_output):
        diagVar.append(np.zeros(count,dtype=dtype,order='F'))
    for vi,vn in enumerate(diag_output):
        bpFileReader.Get(diag_get[vi], diagVar[vi], adios2.Mode.Sync)
    bpFileReader.EndStep()

    #I_out data
    print("start Iout")
    bpFileReader.BeginStep()
    #if xy_rank>11:
    #if yz_procs==nprocs or (yz_procs!=nprocs and yz_rank ==0):
    #if worker_comm.rank ==3 or worker_comm.rank==2 :
    var_get16 = bpIO.InquireVariable(f"I_out.{begin_iter}")
    if var_get16 is None:
        print(comm.rank,"var_get16 is not working")
    var_get16.SetSelection([start4,count4])
    inVar16=np.zeros(count4,dtype=dtype,order='F')
    bpFileReader.Get(var_get16, inVar16, adios2.Mode.Sync)
    bpFileReader.EndStep()
    #if yz_procs==nprocs or (yz_procs!=nprocs and yz_rank ==0):
    #if worker_comm.rank ==3 or worker_comm.rank==2 :
    #if xy_rank>11:
    if multifile:
        prim2_I.I_out[start4[0]:end4[0],start4[1]:end4[1]]=inVar16
    else:
        prim2_I.I_out[i,start4[0]:end4[0],start4[1]:end4[1]]=inVar16
    worker_comm.Barrier()

    #tau_slice data
    #if yz_procs==nprocs or (yz_procs!=nprocs and yz_rank ==0):
    #if worker_comm.rank ==0 or worker_comm.rank==1 :
    inVar12=[]
    for c,li in enumerate(tau_lev):
        for vi in range(tau_nslvar):
            ind=c*tau_nslvar+vi
            inVar12.append(np.zeros(count4,dtype=dtype,order='F'))
    bpFileReader.BeginStep()
    if xcol_rank ==0:
        var_get12=[]
        var_time2d_get = bpIO.InquireVariable("Run.time.2d")
        for c,li in enumerate(tau_lev):
            for vi in range(tau_nslvar):
                v_str=f"tau_slice_{li}.{vi}.{begin_iter}"
                #print(v_str)
                var_get12.append( bpIO.InquireVariable(v_str))
        if var_get12[0] is None:
            print(comm.rank,"var_get12 inquire is not working")
        inVarTime2d=np.zeros((1),dtype='f4')
        if var_time2d_get is None:
            print(comm.rank,"var_time2d_get is not working")
        for c,li in enumerate(tau_lev):
            for vi in range(tau_nslvar):
                ind=c*tau_nslvar+vi
                var_get12[ind].SetSelection([start4,count4])
                #inVar12.append(np.zeros(count4,dtype=dtype,order='F'))
        bpFileReader.Get(var_time2d_get, inVarTime2d, adios2.Mode.Sync)
        for c,li in enumerate(tau_lev):
            for vi in range(tau_nslvar):
                ind=c*tau_nslvar+vi
                bpFileReader.Get(var_get12[ind], inVar12[ind], adios2.Mode.Sync)
    err=bpFileReader.EndStep()
    #if yz_procs==nprocs or (yz_procs!=nprocs and yz_rank ==0):
    #if worker_comm.rank ==0 or worker_comm.rank==1 :
    if xcol_rank ==0:
        if multifile:
            for c,li in enumerate(tau_lev):
                for vi in range(tau_nslvar):
                    ind=c*tau_nslvar+vi
                    prim2_tau[f'tau_slice_{li}_{vi}'][start4[0]:end4[0],start4[1]:end4[1]]=inVar12[ind]
            if yz_rank ==0 :
                prim2_tau.attrs['run_time_2d']=inVarTime2d[0]
    print("tau_slice IO done")
    worker_comm.Barrier()
 
    #yz_slice IO
    if yz_slice:
        bpFileReader.BeginStep()
        var_get_yz={}
        for vi in format_yz_lev:
            #if str(rank) in yz_rank_info[f'yz_slice_{vi}'] and rank<=yz_rank_info[vi][1]:
            if str(rank) in rank_info[f'yz_slice_{vi}']:
                var_name=f"yz_slice_{vi}.{begin_iter}"
                print('var_name ',var_name, vi)
                var_get_yz[var_name]=bpIO.InquireVariable(var_name)
        inVar_yz={}
        for c,vi in var_get_yz.items():
            print('var_get_yz ',c,type(vi))
            vi.SetSelection([start_yz,count_yz])
            inVar_yz[c]=np.zeros(count_yz,dtype=dtype,order='C')
        for c,vi in var_get_yz.items():
            bpFileReader.Get(vi, inVar_yz[c], adios2.Mode.Deferred)
            print('inVar_yz ',c,type(inVar_yz[c]))
        bpFileReader.EndStep()

        print("yz_slice IO done")
        for c,vi in inVar_yz.items():
            print(rank,'write inVar_yz ',c,inVar_yz[c][0,0,0])
        if multifile:
            for c,vi in inVar_yz.items():
                prim2_yz[c[:13]][start_yz[0]:end_yz[0],start_yz[1]:end_yz[1],start_yz[2]:end_yz[2]]=inVar_yz[c]

    #xy_slice IO
    bpFileReader.BeginStep()
    var_get_xy={}
    for vi in format_xy_lev:
        if str(rank) in rank_info[f'xy_slice_{vi}']:
            var_name=f"xy_slice_{vi}.{begin_iter}"
            var_get_xy[var_name]=bpIO.InquireVariable(var_name)
    inVar_xy={} 
    for c,vi in var_get_xy.items():
        vi.SetSelection([start3,count3])
        inVar_xy[c]=np.zeros(count3,dtype=dtype,order='C')
    for c,vi in var_get_xy.items():
        bpFileReader.Get(vi, inVar_xy[c], adios2.Mode.Sync)
    bpFileReader.EndStep()
    print("xy_slice IO done")
    if multifile:
        for c,vi in inVar_xy.items():
            prim2_xy[c[:13]][start3[0]:end3[0],start3[1]:end3[1],start3[2]:end3[2]]=inVar_xy[c]
        #if worker_comm.rank ==0 :
        #    prim2_xy.attrs['run_time_2d']=inVarTime2d[0]
    #else:
    #    prim2_xy.xy_slice_0032[i,start3[0]:end3[0],start3[1]:end3[1],start3[2]:end3[2]]=inVar15



    #corona emission variables
    #if bpFileReader.BeginStep()==adios2.StepStatus.OK:
    if dem:
        bpFileReader.BeginStep()
        if ycol_rank == 0:
            #print('corona ycol_rank = 0, rank=',comm.rank)
            var_get21=bpIO.InquireVariable("corona_emission_adj_fil_y.000000")
            var_get22=bpIO.InquireVariable("corona_emission_adj_dem_y.000000")
            var_get23=bpIO.InquireVariable("corona_emission_adj_vlos_y.000000")
            var_get24=bpIO.InquireVariable("corona_emission_adj_vrms_y.000000")
            var_get21.SetSelection([start2[0],count2[0]])
            print(rank," after setselection coro 22")
            var_get22.SetSelection([start2[0],count2[0]])
            var_get23.SetSelection([start2[0],count2[0]])
            var_get24.SetSelection([start2[0],count2[0]])
            inVar21=np.zeros(count2[0],dtype=dtype,order='C')
            inVar22=np.zeros(count2[0],dtype=dtype,order='C')
            inVar23=np.zeros(count2[0],dtype=dtype,order='C')
            inVar24=np.zeros(count2[0],dtype=dtype,order='C')
            bpFileReader.Get(var_get21, inVar21,adios2.Mode.Sync)
            bpFileReader.Get(var_get22, inVar22,adios2.Mode.Sync)
            bpFileReader.Get(var_get23, inVar23,adios2.Mode.Sync)
            bpFileReader.Get(var_get24, inVar24,adios2.Mode.Sync)
            '''
            if ycol_rank == 0:
                if multifile:
                    prim2_corona.corona_y1[start2[0][0]:end2[0][0],start2[0][1]:end2[0][1],start2[0][2]:end2[0][2]]=inVar21
                    prim2_corona.corona_y2[start2[0][0]:end2[0][0],start2[0][1]:end2[0][1],start2[0][2]:end2[0][2]]=inVar22
                    prim2_corona.corona_y3[start2[0][0]:end2[0][0],start2[0][1]:end2[0][1],start2[0][2]:end2[0][2]]=inVar23
                    prim2_corona.corona_y4[start2[0][0]:end2[0][0],start2[0][1]:end2[0][1],start2[0][2]:end2[0][2]]=inVar24
            '''

        bpFileReader.EndStep()
        bpFileReader.BeginStep()
        if xcol_rank == 0:
            print('corona xcol_rank = 0, rank=',comm.rank)
            var_get17=bpIO.InquireVariable("corona_emission_adj_fil_x.000000")
            var_get18=bpIO.InquireVariable("corona_emission_adj_dem_x.000000")
            var_get19=bpIO.InquireVariable("corona_emission_adj_vlos_x.000000")
            var_get20=bpIO.InquireVariable("corona_emission_adj_vrms_x.000000")
            #print("before setselection coro 17",start2[0],start2[1],start2[2],count2[0],count2[1],count2[2],comm.rank)
            var_get17.SetSelection([start2[1],count2[1]])
            print("after setselection coro 17")
            var_get18.SetSelection([start2[1],count2[1]])
            var_get19.SetSelection([start2[1],count2[1]])
            var_get20.SetSelection([start2[1],count2[1]])
            inVar17=np.zeros(count2[1],dtype=dtype,order='C')
            inVar18=np.zeros(count2[1],dtype=dtype,order='C')
            inVar19=np.zeros(count2[1],dtype=dtype,order='C')
            inVar20=np.zeros(count2[1],dtype=dtype,order='C')
            bpFileReader.Get(var_get17, inVar17,adios2.Mode.Sync)
            bpFileReader.Get(var_get18, inVar18,adios2.Mode.Sync)
            bpFileReader.Get(var_get19, inVar19,adios2.Mode.Sync)
            bpFileReader.Get(var_get20, inVar20,adios2.Mode.Sync)
        '''
        if xcol_rank == 0:
            if multifile:
                prim2_corona.corona_x1[start2[1][0]:end2[1][0],start2[1][1]:end2[1][1],start2[1][2]:end2[1][2]]=inVar17
                prim2_corona.corona_x2[start2[1][0]:end2[1][0],start2[1][1]:end2[1][1],start2[1][2]:end2[1][2]]=inVar18
                prim2_corona.corona_x3[start2[1][0]:end2[1][0],start2[1][1]:end2[1][1],start2[1][2]:end2[1][2]]=inVar19
                prim2_corona.corona_x4[start2[1][0]:end2[1][0],start2[1][1]:end2[1][1],start2[1][2]:end2[1][2]]=inVar20
         
        '''
        bpFileReader.EndStep()
        bpFileReader.BeginStep()
        if zcol_rank == 0:
            print("before inquire coro 25")
            var_get25=bpIO.InquireVariable("corona_emission_adj_fil_z.000000")
            var_get26=bpIO.InquireVariable("corona_emission_adj_dem_z.000000")
            var_get27=bpIO.InquireVariable("corona_emission_adj_vlos_z.000000")
            var_get28=bpIO.InquireVariable("corona_emission_adj_vrms_z.000000")
            #print("before setselection coro 25",start2[0],start2[1],start2[2],count2[0],count2[1],count2[2],comm.rank)
            var_get25.SetSelection([start2[2],count2[2]])
            print("after setselection coro 25")
            var_get26.SetSelection([start2[2],count2[2]])
            var_get27.SetSelection([start2[2],count2[2]])
            var_get28.SetSelection([start2[2],count2[2]])
            inVar25=np.zeros(count2[2],dtype=dtype,order='C')
            inVar26=np.zeros(count2[2],dtype=dtype,order='C')
            inVar27=np.zeros(count2[2],dtype=dtype,order='C')
            inVar28=np.zeros(count2[2],dtype=dtype,order='C')
            bpFileReader.Get(var_get25, inVar25,adios2.Mode.Sync)
            bpFileReader.Get(var_get26, inVar26,adios2.Mode.Sync)
            bpFileReader.Get(var_get27, inVar27,adios2.Mode.Sync)
            bpFileReader.Get(var_get28, inVar28,adios2.Mode.Sync)
        bpFileReader.EndStep()
        '''
        if zcol_rank == 0:
            if multifile:
                prim2_corona.corona_z1[start2[2][0]:end2[2][0],start2[2][1]:end2[2][1],start2[2][2]:end2[2][2]]=inVar25
                prim2_corona.corona_z2[start2[2][0]:end2[2][0],start2[2][1]:end2[2][1],start2[2][2]:end2[2][2]]=inVar26
                prim2_corona.corona_z3[start2[2][0]:end2[2][0],start2[2][1]:end2[2][1],start2[2][2]:end2[2][2]]=inVar27
                prim2_corona.corona_z4[start2[2][0]:end2[2][0],start2[2][1]:end2[2][1],start2[2][2]:end2[2][2]]=inVar28
        '''

    worker_comm.Barrier()
    if multifile:
        start_time=MPI.Wtime()
        #print("prim_0 write to file")
        prim2.prim_0[start[0]:end[0],start[1]:end[1],start[2]:end[2]]=inVar
        prim2.prim_1[start[0]:end[0],start[1]:end[1],start[2]:end[2]]=inVar1
        prim2.prim_2[start[0]:end[0],start[1]:end[1],start[2]:end[2]]=inVar2
        prim2.prim_3[start[0]:end[0],start[1]:end[1],start[2]:end[2]]=inVar3
        prim2.prim_4[start[0]:end[0],start[1]:end[1],start[2]:end[2]]=inVar4
        prim2.prim_5[start[0]:end[0],start[1]:end[1],start[2]:end[2]]=inVar5
        prim2.prim_6[start[0]:end[0],start[1]:end[1],start[2]:end[2]]=inVar6
        prim2.prim_7[start[0]:end[0],start[1]:end[1],start[2]:end[2]]=inVar7
        prim2.prim_8[start[0]:end[0],start[1]:end[1],start[2]:end[2]]=inVar8
    for vi,vn in enumerate(eos_output):
        prim2[vn][start[0]:end[0],start[1]:end[1],start[2]:end[2]]=eosVar[vi]
    for vi,vn in enumerate(diag_output):
        prim2[vn][start[0]:end[0],start[1]:end[1],start[2]:end[2]]=diagVar[vi]
    print('Total 3D write time:',comm.rank,MPI.Wtime()-start_time)

    if need_zip:
        store.close()
        store_I.close()
        store_xy.close()
        store_yz.close()
        store_tau.close()
        store_corona.close()
        

bpFileReader.Close()

