#!/bin/bash

module load cray-python/3.10.10

start_idx=$(( $NUM_NODES_PER_BATCH * ($CURRENT_BATCH - 1) ))
end_idx=$(( $NUM_NODES_PER_BATCH * $CURRENT_BATCH ))

srun -n $NUM_NODES_PER_BATCH -c 16 python $PY_FILE -s $start_idx -e $end_idx
'


