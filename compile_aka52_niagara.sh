module purge
module load NiaEnv/2019b
module load gcc/9.4.0
module load openmpi/4.1.1
module load hdf5/1.10.9
module load python/2.7.15


export LD_LIBRARY_PATH=/scinet/niagara/software/2019b/opt/base/python/2.7.15/lib:$LD_LIBRARY_PATH
export PYTHONPATH=/scinet/niagara/software/2019b/opt/base/python/2.7.15/bin/python2.7

#export OMP_NUM_THREADS=1
#export OMP_SCHEDULE=dynamic
#export OMP_PROC_BIND=true
#export OMPI_MCA_btl_portals4_use_rdma=0
#
## For MPI-tags:
#export MPIR_CVAR_CH4_OFI_TAG_BITS=26
#export MPIR_CVAR_CH4_OFI_RANK_BITS=13


