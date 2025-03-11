> TODO

> not sure where to put below (snippet from DiRAC)

gpu_launch.sh wrapper script
The gpu_launch.sh wrapper script is required to set the correct binding of GPU to MPI processes and the correct binding of interconnect interfaces to MPI process and GPU. We provide this centrally for convenience but its contents are simple:

```
#!/bin/bash

# Compute the raw process ID for binding to GPU and NIC
lrank=$((SLURM_PROCID % SLURM_NTASKS_PER_NODE))

# Bind the process to the correct GPU and NIC
export CUDA_VISIBLE_DEVICES=${lrank}
export UCX_NET_DEVICES=mlx5_${lrank}:1
```