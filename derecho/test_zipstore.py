import zarr
import xarray as xr
from os import system
import time

#mkdir 7z
#mv 7z2407-linux-x64.tar.xz 7z/
#cd 7z
#tar -xf 7z2407-linux-x64.tar.xz

cmds="/glade/work/haiyingx/7z/7zz  a -tzip I_out.2.lossless.single.122.zarr.zip I_out.2.lossless.single.122.zarr/."
system(cmds)
filename='I_out.2.lossless.single.122.zarr.zip'
#file=zarr.group(store)
start = time.perf_counter()
store=zarr.ZipStore(filename)
print(xr.open_zarr(store))
print("zipstore open time = ",time.perf_counter()-start)
start = time.perf_counter()
print(xr.open_zarr('I_out.2.lossless.single.122.zarr'))
print("direct zarr open time = ",time.perf_counter()-start)
