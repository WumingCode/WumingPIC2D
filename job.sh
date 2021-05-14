#!/bin/sh
#PJM --name "em2d_test"
#PJM -L  "node=1"                      # 割当ノード数 12 （3次元形状）
#PJM -L  "rscgrp=small"                    # リソースグループの指定
#PJM -L  "elapse=01:00:00"                 # 経過時間制限 1時間
#PJM --mpi "max-proc-per-node=4"           # 1ノードあたりに生成するMPIプロセス数の上限値
#PJM -o sys_out
#PJM -e sys_err

export OMP_NUM_THREADS=12        # スレッド数の指定
export PLE_MPI_STD_EMPTYFILE=OFF

. /vol0004/apps/oss/spack/share/spack/setup-env.sh
spack load /7cxah5f
export LD_LIBRARY_PATH=/lib64:$LD_LIBRARY_PATH

mpiexec -n 4 -stdout sys_out -stderr sys_err ./em2dsk_test.out

