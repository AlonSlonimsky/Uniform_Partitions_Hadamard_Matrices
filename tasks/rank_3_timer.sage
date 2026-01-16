load("../collections/class_representatives.sage")
load("../tools/generating_matrices.sage")
load("../tools/are_hadamard_equivalent.sage")
load("../tools/hadamard_tools.sage")

import multiprocessing
import time
import sys
from datetime import datetime

# This file provides a sample for doing searches for specific partitions
# Can change which order of hadamard matrix by changing the first line in long_job
# Can change which algorithm is used by changing the function
# Can change the timeout and pool settings in the bottom 


def long_job(i):
    with open(f"file{i}.txt", "a") as f:
        mat = classes24[i] # which order of H(n)
        
        f.write("------\n")
        print("------")
        
        f.write(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + "\n")
        print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        
        f.write(f"Matrix: {i}\n")
        print(f"Matrix: {i}")
        
        f.write(f"{mat}\n")
        print(mat)
        
        f.write("Partitioned:\n")
        print("Partitioned:")

        f.flush() # force the write to the file to be completed now
        
        result = rank_3_partition_24(mat) # which algorithm/function

        if result == None: # No found partition
            print("Completed Search")
            return

        result = partition_3x3(result)

        print(are_hadamard_equivalent(mat, result)) # check just in case

        print(get_matrix_of_3x3_ranks(result)) # another check

        signature = get_matrix_signature(result)

        f.write(f"{result}\n")
        print(result)

        f.write("Signature:\n")
        print("Signature")
        f.write(f"{signature}\n")
        print(signature)
        
        f.write(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + "\n") # completion time with formatting
        print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    

# necessary for multiprocessing
# effectively just identifies the main process
if __name__ == "__main__":     
    inputs = range(60)

    timeout_seconds = 60*60*24*2  # 2 full days
    
    # Argument in Pool() is how many processes to run at once. Leave empty for default, which should be total number of cores.
    with multiprocessing.Pool(1) as pool: 
        result = pool.map_async(long_job, inputs)
        try:
            result.get(timeout=timeout_seconds)
        except multiprocessing.TimeoutError:
            pool.terminate()  # terminate all running tasks once timeout seconds passed
        finally:
            pool.join()
    