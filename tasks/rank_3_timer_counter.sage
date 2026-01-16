load("../collections/class_representatives.sage")
load("../tools/generating_matrices.sage")
load("../tools/are_hadamard_equivalent.sage")
load("../tools/hadamard_tools.sage")

import multiprocessing
import time
import sys
from datetime import datetime

# I *think* that this code will work on the distributed system.
# It starts a bunch of processes then just waits 2 days before force terminating the remaining ones
# With enough cores most should complete quickly, but there may be some stubborn ones


def long_job(arg_tuple):
    row2 = arg_tuple[0]
    row3 = arg_tuple[1]
    mat = arg_tuple[2]

    count = rank_3_partition_24_count(mat, row2, row3)
    with open(f"{row2}-{row3}.txt", "a") as f:
        f.write(f"Count: {count}\n")
    

# necessary for multiprocessing
# effectively just identifies the main process
if __name__ == "__main__":     

    mat = classes24[36] # 36 is the mostly arbitrary class we've chosen for a complete search

    inputs = [(row2, row3, mat) for row2, row3 in itertools.combinations(range(1, 24), 2)]

    timeout_seconds = 60*60*24*2  # 2 full days
        
    with multiprocessing.Pool(1) as pool: # NOTE: will melt pc if run over all cores
        result = pool.map_async(long_job, inputs)
        try:
            result.get(timeout=timeout_seconds)
        except multiprocessing.TimeoutError:
            pool.terminate()  # terminate all running tasks once timeout seconds passed
        finally:
            pool.join()
    