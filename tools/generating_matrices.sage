load("../tools/hadamard_tools.sage")
import random
import itertools
import copy

#------------------------------------------------------
# random_permutation_matrix
#
# PURPOSE: Generate a random nxn permutation matrix.
# PARAMETERS:
#     int n: The order of the permutation matrix.
# Returns: The random nxn permuation matrix.
#------------------------------------------------------
def random_permutation_matrix(n: int):
    permutation = list(range(n))
    random.shuffle(permutation)
    
    mat = matrix(n, n, 0) # 0 matrix

    for i in range(n):
        mat[i, permutation[i]] = 1
    
    return mat

#------------------------------------------------------
# rank_n_random_partition
#
# PURPOSE: Generate (or attempt to) a random rank n 3x3 partition of the given matrix by performing random permutations till a valid one is found.
#          Not deterministic, may never stop.
# PARAMETERS:
#     Matrix mat: The matrix to be permuted.
#     int n: The desired rank of the result.
# Returns: If this function returns, it returns the permutation of mat in a rank n 3x3 partition.
#------------------------------------------------------
def rank_n_random_partition(mat: Matrix, n: int) -> Matrix:
    while True: # forever

        #NOTE: technically the possible permuation matrices could be restricted, as some permutations are redundant,
        #      i.e. swapping the last two rows will not affect the 3x3 partition's rank.
        #      This is not however a very impressive improvement, and does not resolve the O(n!) scaling.
        first = random_permutation_matrix(mat.nrows())
        second = random_permutation_matrix(mat.nrows())

        equivalent_mat = first * mat * second # permute mat

        sub_blocks = get_3x3_sub_blocks(equivalent_mat)

        bad_sub_block_found = False
        for sub_block in sub_blocks:
            if sub_block.rank() != n:
                bad_sub_block_found = True
                break
        if bad_sub_block_found:
            continue # try next random pair of permutation matrices

        return partition_3x3(equivalent_mat) # parititon the new matrix and return it
    
    return None # unreachable


#------------------------------------------------------
# rank_3_partition_24
#
# PURPOSE: Generate a rank 3 partition of an H(24) while restricting the search size by enforcing the edges to be of a valid type
# PARAMETERS:
#     Matrix mat: The matrix to be permuted.
# Returns: The rank 3 partition, or None if there isn't one
#------------------------------------------------------
def rank_3_partition_24(A: Matrix) -> Matrix:
    # avoid modifying A directly
    B = copy.deepcopy(A)

    # verify matrix is normalized along the first row
    for c in range(A.ncols()):
        if B[0,c] == -1:
            B[:,c] = -B[:,c] # negate the entire column

    # decide row 2 and 3, row 1 is fixed
    for row2, row3 in itertools.combinations(range(1,24), 2):

        C = copy.deepcopy(B)
        # these swaps allows all row decisions to be from 3-24
        C.swap_rows(1,row2)
        C.swap_rows(2,row3)

        A_cols_1 = []
        B_cols_1 = []
        C_cols_1 = []
        D_cols_1 = []

        # get the type of every column from 1-24
        for i in range(24):
            element_triple = (C[0,i],C[1,i],C[2,i])
            if element_triple == (1,1,1) or element_triple == (-1,-1,-1):
                A_cols_1.append(i)
            elif element_triple == (-1,1,1) or element_triple == (1,-1,-1):
                B_cols_1.append(i)
            elif element_triple == (1,-1,1) or element_triple == (-1,1,-1):
                C_cols_1.append(i)
            elif element_triple == (1,1,-1) or element_triple == (-1,-1,1):
                D_cols_1.append(i)

        # fix the a column positions

        a_col_1 = A_cols_1[0]
        a_col_2 = A_cols_1[1]
        a_col_3 = A_cols_1[2]
        a_col_4 = A_cols_1[3]
        a_col_5 = A_cols_1[4]
        a_col_6 = A_cols_1[5]

        # these two arrays and their updated versions keep track of our current partial rank 3 partition

        decided_rows_2 = [0,1,2] # fixed
        decided_cols_1 = []

        # decide col 1,2,3 (a,b,c)
        for b_col_1, c_col_1 in itertools.product(B_cols_1, C_cols_1):
            B_cols_2 = B_cols_1.copy(); B_cols_2.remove(b_col_1)
            C_cols_2 = C_cols_1.copy(); C_cols_2.remove(c_col_1)

            decided_cols_2 = decided_cols_1 + [a_col_1, b_col_1, c_col_1]

            # A,B,C starts from 2 here as the first 3 are fixed
            A_rows_2 = []
            B_rows_2 = []
            C_rows_2 = []
            D_rows_1 = []

            # categorize rows as we did the columns
            # skip the first three rows, they are forces to be of types a,b,c
            for i in range(3,24):
                element_triple = (C[i,a_col_1],C[i,b_col_1],C[i,c_col_1])
                if element_triple == (1,1,1) or element_triple == (-1,-1,-1):
                    A_rows_2.append(i)
                elif element_triple == (-1,1,1) or element_triple == (1,-1,-1):
                    B_rows_2.append(i)
                elif element_triple == (1,-1,1) or element_triple == (-1,1,-1):
                    C_rows_2.append(i)
                elif element_triple == (1,1,-1) or element_triple == (-1,-1,1):
                    D_rows_1.append(i)

            # fix the order of the a rows
            a_row_2 = A_rows_2[0]
            a_row_3 = A_rows_2[1]
            a_row_4 = A_rows_2[2]
            a_row_5 = A_rows_2[3]
            a_row_6 = A_rows_2[4]


            for b_col_2, c_col_2 in itertools.product(B_cols_2, C_cols_2):
                B_cols_3 = B_cols_2.copy(); B_cols_3.remove(b_col_2)
                C_cols_3 = C_cols_2.copy(); C_cols_3.remove(c_col_2)

                decided_cols_3 = decided_cols_2 + [a_col_2, b_col_2, c_col_2]

                # this choice is always valid, so no check & continue

                for b_row_2, c_row_2 in itertools.product(B_rows_2, C_rows_2):
                    B_rows_3 = B_rows_2.copy(); B_rows_3.remove(b_row_2)
                    C_rows_3 = C_rows_2.copy(); C_rows_3.remove(c_row_2)

                    decided_rows_3 = decided_rows_2 + [a_row_2, b_row_2, c_row_2]

                    if not check_row(C, a_row_2, b_row_2, c_row_2, decided_cols_3, 2):
                        continue # if this wasn't valid then continue
                        # continue goes to the next iteration of the above loop

                    for b_col_3, d_col_1 in itertools.product(B_cols_3, D_cols_1):
                        B_cols_4 = B_cols_3.copy(); B_cols_4.remove(b_col_3)
                        D_cols_2 = D_cols_1.copy(); D_cols_2.remove(d_col_1)
                        
                        decided_cols_4 = decided_cols_3 + [a_col_3, b_col_3, d_col_1]

                        if not check_column(C, a_col_3, b_col_3, d_col_1, decided_rows_3, 2):
                            continue

                        for b_row_3, d_row_1 in itertools.product(B_rows_3, D_rows_1):
                            B_rows_4 = B_rows_3.copy(); B_rows_4.remove(b_row_3)
                            D_rows_2 = D_rows_1.copy(); D_rows_2.remove(d_row_1)
                            
                            decided_rows_4 = decided_rows_3 + [a_row_3, b_row_3, d_row_1]

                            if not check_row(C, a_row_3, b_row_3, d_row_1, decided_cols_4, 3):
                                continue

                            for b_col_4, d_col_2 in itertools.product(B_cols_4, D_cols_2):
                                B_cols_5 = B_cols_4.copy(); B_cols_5.remove(b_col_4)
                                D_cols_3 = D_cols_2.copy(); D_cols_3.remove(d_col_2)

                                decided_cols_5 = decided_cols_4 + [a_col_4, b_col_4, d_col_2]

                                if not check_column(C, a_col_4, b_col_4, d_col_2, decided_rows_4, 3):
                                    continue
                                
                                # the remaining b's can be fixed here to avoid duplication further on
                                b_col_5 = B_cols_5.pop(0)
                                b_col_6 = B_cols_5.pop(0)

                                for b_row_4, d_row_2 in itertools.product(B_rows_4, D_rows_2):
                                    B_rows_5 = B_rows_4.copy(); B_rows_5.remove(b_row_4)
                                    D_rows_3 = D_rows_2.copy(); D_rows_3.remove(d_row_2)
                                    
                                    decided_rows_5 = decided_rows_4 + [a_row_4, b_row_4, d_row_2]

                                    if not check_row(C, a_row_4, b_row_4, d_row_2, decided_cols_5, 4):
                                        continue

                                    # the remaining b's can be fixed here to avoid duplication further on
                                    b_row_5 = B_rows_5.pop(0)
                                    b_row_6 = B_rows_5.pop(0)

                                    for c_col_3, d_col_3 in itertools.product(C_cols_3, D_cols_3):
                                        C_cols_4 = C_cols_3.copy(); C_cols_4.remove(c_col_3)
                                        D_cols_4 = D_cols_3.copy(); D_cols_4.remove(d_col_3)

                                        decided_cols_6 = decided_cols_5 + [a_col_5, c_col_3, d_col_3]

                                        if not check_column(C, a_col_5, c_col_3, d_col_3, decided_rows_5, 4):
                                            continue

                                        for c_row_3, d_row_3 in itertools.product(C_rows_3, D_rows_3):
                                            C_rows_4 = C_rows_3.copy(); C_rows_4.remove(c_row_3)
                                            D_rows_4 = D_rows_3.copy(); D_rows_4.remove(d_row_3)

                                            decided_rows_6 = decided_rows_5 + [a_row_5, c_row_3, d_row_3]

                                            if not check_row(C, a_row_5, c_row_3, d_row_3, decided_cols_6, 5):
                                                continue

                                            for c_col_4, d_col_4 in itertools.product(C_cols_4, D_cols_4):
                                                C_cols_5 = C_cols_4.copy(); C_cols_5.remove(c_col_4)
                                                D_cols_5 = D_cols_4.copy(); D_cols_5.remove(d_col_4)

                                                decided_cols_7 = decided_cols_6 + [a_col_6, c_col_4, d_col_4]

                                                if not check_column(C, a_col_6, c_col_4, d_col_4, decided_rows_6, 5):
                                                    continue

                                                for c_row_4, d_row_4 in itertools.product(C_rows_4, D_rows_4):
                                                    C_rows_5 = C_rows_4.copy(); C_rows_5.remove(c_row_4)
                                                    D_rows_5 = D_rows_4.copy(); D_rows_5.remove(d_row_4)

                                                    decided_rows_7 = decided_rows_6 + [a_row_6, c_row_4, d_row_4]

                                                    if not check_row(C, a_row_6, c_row_4, d_row_4, decided_cols_7, 6):
                                                        continue

                                                    for c_col_5, d_col_5 in itertools.product(C_cols_5, D_cols_5):
                                                        C_cols_6 = C_cols_5.copy(); C_cols_6.remove(c_col_5)
                                                        D_cols_6 = D_cols_5.copy(); D_cols_6.remove(d_col_5)

                                                        decided_cols_8 = decided_cols_7 + [b_col_5, c_col_5, d_col_5]

                                                        if not check_column(C, b_col_5, c_col_5, d_col_5, decided_rows_7, 6):
                                                            continue

                                                        for c_row_5, d_row_5 in itertools.product(C_rows_5, D_rows_5):
                                                            C_rows_6 = C_rows_5.copy(); C_rows_6.remove(c_row_5)
                                                            D_rows_6 = D_rows_5.copy(); D_rows_6.remove(d_row_5)

                                                            decided_rows_8 = decided_rows_7 + [b_row_5, c_row_5, d_row_5]

                                                            if not check_row(C, b_row_5, c_row_5, d_row_5, decided_cols_8, 7):
                                                                continue

                                                            for c_col_6, d_col_6 in itertools.product(C_cols_6, D_cols_6):
                                                                C_cols_7 = C_cols_6.copy(); C_cols_7.remove(c_col_6)
                                                                D_cols_7 = D_cols_6.copy(); D_cols_7.remove(d_col_6)

                                                                decided_cols_9 = decided_cols_8 + [b_col_6, c_col_6, d_col_6]

                                                                if not check_column(C, b_col_6, c_col_6, d_col_6, decided_rows_8, 7):
                                                                    continue

                                                                for c_row_6, d_row_6 in itertools.product(C_rows_6, D_rows_6):
                                                                    C_rows_7 = C_rows_6.copy(); C_rows_7.remove(c_row_6)
                                                                    D_rows_7 = D_rows_6.copy(); D_rows_7.remove(d_row_6)

                                                                    decided_rows_9 = decided_rows_8 + [b_row_6, c_row_6, d_row_6]

                                                                    if not check_row(C, b_row_6, c_row_6, d_row_6, decided_cols_9, 8):
                                                                        continue
                                                                    
                                                                    # else, Good, YAY!
                                                                    # perform the transformation in place
                                                                    return matrix(C.nrows(), C.ncols(), lambda i, j: C[decided_rows_9[i], decided_cols_9[j]])
    return None

#------------------------------------------------------
# rank_3_partition_24
#
# PURPOSE: Count the number of rank 3 partitions our algorithm finds. The choice of second and third row made initially is made an argument
# PARAMETERS:
#     Matrix mat: The matrix to be permuted.
#     int row2: The second row.
#     int row3: The third row.
# Returns: The number of found rank 3 partitions
#------------------------------------------------------
def rank_3_partition_24_count(A: Matrix, row2: int, row3: int) -> Matrix:
    count = 0

    # avoid modifying A directly
    B = copy.deepcopy(A)

    # verify matrix is normalized along the first row
    for c in range(A.ncols()):
        if B[0,c] == -1:
            B[:,c] = -B[:,c] # negate the entire column

    # decide row 2 and 3, row 1 is fixed
    

    C = copy.deepcopy(B)
    # these swaps allows all row decisions to be from 3-24
    C.swap_rows(1,row2)
    C.swap_rows(2,row3)

    A_cols_1 = []
    B_cols_1 = []
    C_cols_1 = []
    D_cols_1 = []

    # get the type of every column from 1-24
    for i in range(24):
        element_triple = (C[0,i],C[1,i],C[2,i])
        if element_triple == (1,1,1) or element_triple == (-1,-1,-1):
            A_cols_1.append(i)
        elif element_triple == (-1,1,1) or element_triple == (1,-1,-1):
            B_cols_1.append(i)
        elif element_triple == (1,-1,1) or element_triple == (-1,1,-1):
            C_cols_1.append(i)
        elif element_triple == (1,1,-1) or element_triple == (-1,-1,1):
            D_cols_1.append(i)

    # fix the a column positions

    a_col_1 = A_cols_1[0]
    a_col_2 = A_cols_1[1]
    a_col_3 = A_cols_1[2]
    a_col_4 = A_cols_1[3]
    a_col_5 = A_cols_1[4]
    a_col_6 = A_cols_1[5]

    # these two arrays and their updated versions keep track of our current partial rank 3 partition

    decided_rows_2 = [0,1,2] # fixed
    decided_cols_1 = []

    # decide col 1,2,3 (a,b,c)
    for b_col_1, c_col_1 in itertools.product(B_cols_1, C_cols_1):
        B_cols_2 = B_cols_1.copy(); B_cols_2.remove(b_col_1)
        C_cols_2 = C_cols_1.copy(); C_cols_2.remove(c_col_1)

        decided_cols_2 = decided_cols_1 + [a_col_1, b_col_1, c_col_1]

        # A,B,C starts from 2 here as the first 3 are fixed
        A_rows_2 = []
        B_rows_2 = []
        C_rows_2 = []
        D_rows_1 = []

        # categorize rows as we did the columns
        # skip the first three rows, they are forces to be of types a,b,c
        for i in range(3,24):
            element_triple = (C[i,a_col_1],C[i,b_col_1],C[i,c_col_1])
            if element_triple == (1,1,1) or element_triple == (-1,-1,-1):
                A_rows_2.append(i)
            elif element_triple == (-1,1,1) or element_triple == (1,-1,-1):
                B_rows_2.append(i)
            elif element_triple == (1,-1,1) or element_triple == (-1,1,-1):
                C_rows_2.append(i)
            elif element_triple == (1,1,-1) or element_triple == (-1,-1,1):
                D_rows_1.append(i)

        # fix the order of the a rows
        a_row_2 = A_rows_2[0]
        a_row_3 = A_rows_2[1]
        a_row_4 = A_rows_2[2]
        a_row_5 = A_rows_2[3]
        a_row_6 = A_rows_2[4]


        for b_col_2, c_col_2 in itertools.product(B_cols_2, C_cols_2):
            B_cols_3 = B_cols_2.copy(); B_cols_3.remove(b_col_2)
            C_cols_3 = C_cols_2.copy(); C_cols_3.remove(c_col_2)

            decided_cols_3 = decided_cols_2 + [a_col_2, b_col_2, c_col_2]

            # this choice is always valid, so no check & continue

            for b_row_2, c_row_2 in itertools.product(B_rows_2, C_rows_2):
                B_rows_3 = B_rows_2.copy(); B_rows_3.remove(b_row_2)
                C_rows_3 = C_rows_2.copy(); C_rows_3.remove(c_row_2)

                decided_rows_3 = decided_rows_2 + [a_row_2, b_row_2, c_row_2]

                if not check_row(C, a_row_2, b_row_2, c_row_2, decided_cols_3, 2):
                    continue # if this wasn't valid then continue
                    # continue goes to the next iteration of the above loop

                for b_col_3, d_col_1 in itertools.product(B_cols_3, D_cols_1):
                    B_cols_4 = B_cols_3.copy(); B_cols_4.remove(b_col_3)
                    D_cols_2 = D_cols_1.copy(); D_cols_2.remove(d_col_1)
                    
                    decided_cols_4 = decided_cols_3 + [a_col_3, b_col_3, d_col_1]

                    if not check_column(C, a_col_3, b_col_3, d_col_1, decided_rows_3, 2):
                        continue

                    for b_row_3, d_row_1 in itertools.product(B_rows_3, D_rows_1):
                        B_rows_4 = B_rows_3.copy(); B_rows_4.remove(b_row_3)
                        D_rows_2 = D_rows_1.copy(); D_rows_2.remove(d_row_1)
                        
                        decided_rows_4 = decided_rows_3 + [a_row_3, b_row_3, d_row_1]

                        if not check_row(C, a_row_3, b_row_3, d_row_1, decided_cols_4, 3):
                            continue

                        for b_col_4, d_col_2 in itertools.product(B_cols_4, D_cols_2):
                            B_cols_5 = B_cols_4.copy(); B_cols_5.remove(b_col_4)
                            D_cols_3 = D_cols_2.copy(); D_cols_3.remove(d_col_2)

                            decided_cols_5 = decided_cols_4 + [a_col_4, b_col_4, d_col_2]

                            if not check_column(C, a_col_4, b_col_4, d_col_2, decided_rows_4, 3):
                                continue
                            
                            # the remaining b's can be fixed here to avoid duplication further on
                            b_col_5 = B_cols_5.pop(0)
                            b_col_6 = B_cols_5.pop(0)

                            for b_row_4, d_row_2 in itertools.product(B_rows_4, D_rows_2):
                                B_rows_5 = B_rows_4.copy(); B_rows_5.remove(b_row_4)
                                D_rows_3 = D_rows_2.copy(); D_rows_3.remove(d_row_2)
                                
                                decided_rows_5 = decided_rows_4 + [a_row_4, b_row_4, d_row_2]

                                if not check_row(C, a_row_4, b_row_4, d_row_2, decided_cols_5, 4):
                                    continue

                                # the remaining b's can be fixed here to avoid duplication further on
                                b_row_5 = B_rows_5.pop(0)
                                b_row_6 = B_rows_5.pop(0)

                                for c_col_3, d_col_3 in itertools.product(C_cols_3, D_cols_3):
                                    C_cols_4 = C_cols_3.copy(); C_cols_4.remove(c_col_3)
                                    D_cols_4 = D_cols_3.copy(); D_cols_4.remove(d_col_3)

                                    decided_cols_6 = decided_cols_5 + [a_col_5, c_col_3, d_col_3]

                                    if not check_column(C, a_col_5, c_col_3, d_col_3, decided_rows_5, 4):
                                        continue

                                    for c_row_3, d_row_3 in itertools.product(C_rows_3, D_rows_3):
                                        C_rows_4 = C_rows_3.copy(); C_rows_4.remove(c_row_3)
                                        D_rows_4 = D_rows_3.copy(); D_rows_4.remove(d_row_3)

                                        decided_rows_6 = decided_rows_5 + [a_row_5, c_row_3, d_row_3]

                                        if not check_row(C, a_row_5, c_row_3, d_row_3, decided_cols_6, 5):
                                            continue

                                        for c_col_4, d_col_4 in itertools.product(C_cols_4, D_cols_4):
                                            C_cols_5 = C_cols_4.copy(); C_cols_5.remove(c_col_4)
                                            D_cols_5 = D_cols_4.copy(); D_cols_5.remove(d_col_4)

                                            decided_cols_7 = decided_cols_6 + [a_col_6, c_col_4, d_col_4]

                                            if not check_column(C, a_col_6, c_col_4, d_col_4, decided_rows_6, 5):
                                                continue

                                            for c_row_4, d_row_4 in itertools.product(C_rows_4, D_rows_4):
                                                C_rows_5 = C_rows_4.copy(); C_rows_5.remove(c_row_4)
                                                D_rows_5 = D_rows_4.copy(); D_rows_5.remove(d_row_4)

                                                decided_rows_7 = decided_rows_6 + [a_row_6, c_row_4, d_row_4]

                                                if not check_row(C, a_row_6, c_row_4, d_row_4, decided_cols_7, 6):
                                                    continue

                                                for c_col_5, d_col_5 in itertools.product(C_cols_5, D_cols_5):
                                                    C_cols_6 = C_cols_5.copy(); C_cols_6.remove(c_col_5)
                                                    D_cols_6 = D_cols_5.copy(); D_cols_6.remove(d_col_5)

                                                    decided_cols_8 = decided_cols_7 + [b_col_5, c_col_5, d_col_5]

                                                    if not check_column(C, b_col_5, c_col_5, d_col_5, decided_rows_7, 6):
                                                        continue

                                                    for c_row_5, d_row_5 in itertools.product(C_rows_5, D_rows_5):
                                                        C_rows_6 = C_rows_5.copy(); C_rows_6.remove(c_row_5)
                                                        D_rows_6 = D_rows_5.copy(); D_rows_6.remove(d_row_5)

                                                        decided_rows_8 = decided_rows_7 + [b_row_5, c_row_5, d_row_5]

                                                        if not check_row(C, b_row_5, c_row_5, d_row_5, decided_cols_8, 7):
                                                            continue

                                                        for c_col_6, d_col_6 in itertools.product(C_cols_6, D_cols_6):
                                                            C_cols_7 = C_cols_6.copy(); C_cols_7.remove(c_col_6)
                                                            D_cols_7 = D_cols_6.copy(); D_cols_7.remove(d_col_6)

                                                            decided_cols_9 = decided_cols_8 + [b_col_6, c_col_6, d_col_6]

                                                            if not check_column(C, b_col_6, c_col_6, d_col_6, decided_rows_8, 7):
                                                                continue

                                                            for c_row_6, d_row_6 in itertools.product(C_rows_6, D_rows_6):
                                                                C_rows_7 = C_rows_6.copy(); C_rows_7.remove(c_row_6)
                                                                D_rows_7 = D_rows_6.copy(); D_rows_7.remove(d_row_6)

                                                                decided_rows_9 = decided_rows_8 + [b_row_6, c_row_6, d_row_6]

                                                                if not check_row(C, b_row_6, c_row_6, d_row_6, decided_cols_9, 8):
                                                                    continue
                                                                
                                                                # else, Good, YAY!
                                                                # perform the transformation in place
                                                                #return matrix(C.nrows(), C.ncols(), lambda i, j: C[decided_rows_9[i], decided_cols_9[j]])
                                                                count += 1
    return count

    

#------------------------------------------------------
# rank_1_partition_24
#
# PURPOSE: Generate a rank 1 partition of an H(24) while restricting the search size by enforcing the edges to be of a valid type
# PARAMETERS:
#     Matrix mat: The matrix to be permuted.
# Returns: The rank 1 partition, or None if there isn't one
#------------------------------------------------------
def rank_1_partition_24(A: Matrix) -> Matrix:
    # avoid modifying A directly
    B = copy.deepcopy(A) 

    # for every possible triple of rows to be the first three
    for row1, row2, row3 in itertools.combinations(range(24), 3):
        C = copy.deepcopy(B)
        # these swaps allows all further row decisions to be from 3-24
        C.swap_rows(0,row1)
        C.swap_rows(1,row2)
        C.swap_rows(2,row3)

        # must do this here as row 1 is changed
        # normalizing along the first row here to force type a to be (1,1,1) not (-1,-1,-1)
        for c in range(C.ncols()):
            if C[0,c] == -1:
                C[:,c] = -C[:,c]#negate
        
        A_cols = []
        B_cols = []
        C_cols = []
        D_cols = []

        # get column types
        for i in range(24):
            element_triple = (C[0,i],C[1,i],C[2,i])
            if element_triple == (1,1,1) or element_triple == (-1,-1,-1):
                A_cols.append(i)
            elif element_triple == (-1,1,1) or element_triple == (1,-1,-1):
                B_cols.append(i)
            elif element_triple == (1,-1,1) or element_triple == (-1,1,-1):
                C_cols.append(i)
            elif element_triple == (1,1,-1) or element_triple == (-1,-1,1):
                D_cols.append(i)
        
        # choose the two groups of As
        for a_col_1, a_col_2, a_col_3 in itertools.combinations(A_cols, 3):
            A_cols_remaining = A_cols.copy()
            A_cols_remaining.remove(a_col_1)
            A_cols_remaining.remove(a_col_2)
            A_cols_remaining.remove(a_col_3)
            a_col_4 = A_cols_remaining.pop(0)
            a_col_5 = A_cols_remaining.pop(0)
            a_col_6 = A_cols_remaining.pop(0)

            # now collect rows by type
            A_rows = []
            B_rows = []
            C_rows = []
            D_rows = []

            # we know the first 3 are of type a
            for i in range(3,24):
                element_triple = (C[i,a_col_1],C[i,a_col_2],C[i,a_col_3])
                if element_triple == (1,1,1) or element_triple == (-1,-1,-1):
                    A_rows.append(i)
                elif element_triple == (-1,1,1) or element_triple == (1,-1,-1):
                    B_rows.append(i)
                elif element_triple == (1,-1,1) or element_triple == (-1,1,-1):
                    C_rows.append(i)
                elif element_triple == (1,1,-1) or element_triple == (-1,-1,1):
                    D_rows.append(i)

            # first three rows are by contruction type a, so collect the last 3

            a_row_4 = A_rows.pop(0)
            a_row_5 = A_rows.pop(0)
            a_row_6 = A_rows.pop(0)

            decided_cols_1 = [a_col_1, a_col_2, a_col_3, a_col_4, a_col_5, a_col_6] # forced rank 1 here
            decided_rows_1 = [0, 1, 2, a_row_4, a_row_5, a_row_6] # forced rank 1 here

            # check the 1x1 3x3 square
            if not check_row(C, a_row_4, a_row_5, a_row_6, decided_cols_1, len(decided_rows_1)/3, 1):
                continue
        
            for b_col_1, b_col_2, b_col_3 in itertools.combinations(B_cols, 3):
                B_cols_remaining = B_cols.copy()
                B_cols_remaining.remove(b_col_1)
                B_cols_remaining.remove(b_col_2)
                B_cols_remaining.remove(b_col_3)
                b_col_4 = B_cols_remaining.pop(0)
                b_col_5 = B_cols_remaining.pop(0)
                b_col_6 = B_cols_remaining.pop(0)

                decided_cols_2 = decided_cols_1 + [b_col_1, b_col_2, b_col_3, b_col_4, b_col_5, b_col_6]

                # since we are deciding two triples of columns at once, check both at once
                if not check_column(C, b_col_1, b_col_2, b_col_3, decided_rows_1, len(decided_rows_1)/3, 1) \
                or not check_column(C, b_col_4, b_col_5, b_col_6, decided_rows_1, len(decided_rows_1)/3, 1):
                    continue

                for b_row_1, b_row_2, b_row_3 in itertools.combinations(B_rows, 3):
                    B_rows_remaining = B_rows.copy()
                    B_rows_remaining.remove(b_row_1)
                    B_rows_remaining.remove(b_row_2)
                    B_rows_remaining.remove(b_row_3)
                    b_row_4 = B_rows_remaining.pop(0)
                    b_row_5 = B_rows_remaining.pop(0)
                    b_row_6 = B_rows_remaining.pop(0)

                    decided_rows_2 = decided_rows_1 + [b_row_1, b_row_2, b_row_3, b_row_4, b_row_5, b_row_6] 

                    if not check_row(C, b_row_1, b_row_2, b_row_3, decided_cols_2, len(decided_cols_2)/3, 1) \
                    or not check_row(C, b_row_4, b_row_5, b_row_6, decided_cols_2, len(decided_cols_2)/3, 1):
                        continue

                    for c_col_1, c_col_2, c_col_3 in itertools.combinations(C_cols, 3):
                        C_cols_remaining = C_cols.copy()
                        C_cols_remaining.remove(c_col_1)
                        C_cols_remaining.remove(c_col_2)
                        C_cols_remaining.remove(c_col_3)
                        c_col_4 = C_cols_remaining.pop(0)
                        c_col_5 = C_cols_remaining.pop(0)
                        c_col_6 = C_cols_remaining.pop(0)

                        decided_cols_3 = decided_cols_2 + [c_col_1, c_col_2, c_col_3, c_col_4, c_col_5, c_col_6]

                        if not check_column(C, c_col_1, c_col_2, c_col_3, decided_rows_2, len(decided_rows_2)/3, 1) \
                        or not check_column(C, c_col_4, c_col_5, c_col_6, decided_rows_2, len(decided_rows_2)/3, 1):
                            continue

                        for c_row_1, c_row_2, c_row_3 in itertools.combinations(C_rows, 3):
                            C_rows_remaining = C_rows.copy()
                            C_rows_remaining.remove(c_row_1)
                            C_rows_remaining.remove(c_row_2)
                            C_rows_remaining.remove(c_row_3)
                            c_row_4 = C_rows_remaining.pop(0)
                            c_row_5 = C_rows_remaining.pop(0)
                            c_row_6 = C_rows_remaining.pop(0)

                            decided_rows_3 = decided_rows_2 + [c_row_1, c_row_2, c_row_3, c_row_4, c_row_5, c_row_6]

                            if not check_row(C, c_row_1, c_row_2, c_row_3, decided_cols_3, len(decided_cols_3)/3, 1) \
                            or not check_row(C, c_row_4, c_row_5, c_row_6, decided_cols_3, len(decided_cols_3)/3, 1):
                                continue

                            for d_col_1, d_col_2, d_col_3 in itertools.combinations(D_cols, 3):
                                D_cols_remaining = D_cols.copy()
                                D_cols_remaining.remove(d_col_1)
                                D_cols_remaining.remove(d_col_2)
                                D_cols_remaining.remove(d_col_3)
                                d_col_4 = D_cols_remaining.pop(0)
                                d_col_5 = D_cols_remaining.pop(0)
                                d_col_6 = D_cols_remaining.pop(0)

                                decided_cols_4 = decided_cols_3 + [d_col_1, d_col_2, d_col_3, d_col_4, d_col_5, d_col_6]

                                if not check_column(C, d_col_1, d_col_2, d_col_3, decided_rows_3, len(decided_rows_3)/3, 1) \
                                or not check_column(C, d_col_4, d_col_5, d_col_6, decided_rows_3, len(decided_rows_3)/3, 1):
                                    continue

                                for d_row_1, d_row_2, d_row_3 in itertools.combinations(D_rows, 3):
                                    D_rows_remaining = D_rows.copy()
                                    D_rows_remaining.remove(d_row_1)
                                    D_rows_remaining.remove(d_row_2)
                                    D_rows_remaining.remove(d_row_3)
                                    d_row_4 = D_rows_remaining.pop(0)
                                    d_row_5 = D_rows_remaining.pop(0)
                                    d_row_6 = D_rows_remaining.pop(0)

                                    decided_rows_4 = decided_rows_3 + [d_row_1, d_row_2, d_row_3, d_row_4, d_row_5, d_row_6]

                                    if not check_row(C, d_row_1, d_row_2, d_row_3, decided_cols_4, len(decided_cols_4)/3, 1) \
                                    or not check_row(C, d_row_4, d_row_5, d_row_6, decided_cols_4, len(decided_cols_4)/3, 1):
                                        continue
                                    
                                    # partition found!
                                    return matrix(C.nrows(), C.ncols(), lambda i, j: C[decided_rows_4[i], decided_cols_4[j]])
    # didn't find a partition
    return None





#------------------------------------------------------
# check_column
#
# PURPOSE: check if a set of columns have correct rank along a set of rows
# PARAMETERS:
#     Matrix mat: The base matrix
#     int col_1: first column
#     int col_2: second column
#     int col_3: third column
#     list rows: what set of rows to check columns in
#     rank: the desired rank
#     check_edges: Whether we need to check the edges.
#                  For example the rank 1 and 3 algorithms don't need to check the edges because they are rank 1/3 by construction.
# Returns: True if all ranks are good, False otherwise
#------------------------------------------------------
def check_column(mat: Matrix, col_1: int, col_2: int, col_3: int, rows: list, blocks_to_check: int, rank = 3, check_edges = False) -> bool:
    
    # determine whether or not to check the edges
    blocks = range(1, blocks_to_check)
    if check_edges:
        blocks = range(blocks_to_check)

    for i in blocks:
        cur_mat = matrix(3,3,   [mat[rows[i*3],col_1], mat[rows[i*3],col_2], mat[rows[i*3],col_3],
                                mat[rows[i*3+1],col_1], mat[rows[i*3+1],col_2], mat[rows[i*3+1],col_3],
                                mat[rows[i*3+2],col_1], mat[rows[i*3+2],col_2], mat[rows[i*3+2],col_3]])
        if cur_mat.rank() != rank:
            return False

    return True

#------------------------------------------------------
# check_row
#
# PURPOSE: check if a set of rows have correct rank along a set of cols
# PARAMETERS:
#     Matrix mat: The base matrix
#     int row_1: first row
#     int row_2: second row
#     int row_3: third row
#     list cols: what set of cols to check columns in
#     rank: the desired rank
#     check_edges: Whether we need to check the edges.
#                  For example the rank 1 and 3 algorithms don't need to check the edges because they are rank 1/3 by construction.
# Returns: True if all ranks are good, False otherwise
#------------------------------------------------------
def check_row(mat: matrix, row_1: int, row_2: int, row_3: int, cols: list, blocks_to_check: int, rank = 3, check_edges = False) -> bool:

    # determine whether or not to check the edges
    blocks = range(1, blocks_to_check)
    if check_edges:
        blocks = range(blocks_to_check)

    for i in blocks: # starting at 1 would avoid a redundant check
        # technically this would be the transpose, same difference though
        cur_mat = matrix(3,3,   [mat[row_1,cols[i*3]], mat[row_2,cols[i*3]], mat[row_3,cols[i*3]],
                                mat[row_1,cols[i*3+1]], mat[row_2,cols[i*3+1]], mat[row_3,cols[i*3+1]],
                                mat[row_1,cols[i*3+2]], mat[row_2,cols[i*3+2]], mat[row_3,cols[i*3+2]]])
        if cur_mat.rank() != rank:
            return False

    return True

#------------------------------------------------------
# order_20_partition
#
# PURPOSE: Generate a rank n partition of an H(20) while restricting the search size by enforcing the construction to continue only if the previous section is good
# PARAMETERS:
#     Matrix mat: The matrix to be permuted.
#     desired_rank: The desired rank,
# Returns: The matrix, or None if there isn't one
#------------------------------------------------------
def order_20_partition(mat: Matrix, desired_rank=3) -> Matrix:
    cols_available_0 = list(range(20))
    rows_available_0 = list(range(20))

    cols_chosen_1 = []
    rows_chosen_1 = []

    A = mat

    # first col is fixed
    for col0, col1 in itertools.combinations(cols_available_0, 2):
        cols_available_1 = cols_available_0.copy()
        cols_available_1.remove(col0)
        cols_available_1.remove(col1)

        # first row is fixed
        for row0, row1 in itertools.combinations(rows_available_0, 2):
            rows_available_1 = rows_available_0.copy()
            rows_available_1.remove(row0)
            rows_available_1.remove(row1)


            for col2, col3, col4 in itertools.combinations(cols_available_1, 3):
                cols_available_2 = cols_available_1.copy()
                cols_available_2.remove(col2)
                cols_available_2.remove(col3)
                cols_available_2.remove(col4)

                cols_chosen_2 = cols_chosen_1 + [col2, col3, col4]

                for row2, row3, row4 in itertools.combinations(rows_available_1, 3):
                    rows_available_2 = rows_available_1.copy()
                    rows_available_2.remove(row2)
                    rows_available_2.remove(row3)
                    rows_available_2.remove(row4)

                    rows_chosen_2 = rows_chosen_1 + [row2, row3, row4]
                    if not check_row(A, row2, row3, row4, cols_chosen_2, len(cols_chosen_2)/3, rank = desired_rank, check_edges = True):
                        continue

                    for col5, col6, col7 in itertools.combinations(cols_available_2, 3):
                        cols_available_3 = cols_available_2.copy()
                        cols_available_3.remove(col5)
                        cols_available_3.remove(col6)
                        cols_available_3.remove(col7)

                        cols_chosen_3 = cols_chosen_2 + [col5, col6, col7]
                        if not check_column(A, col5, col6, col7, rows_chosen_2, len(rows_chosen_2)/3, rank = desired_rank, check_edges = True):
                            continue

                        for row5, row6, row7 in itertools.combinations(rows_available_2, 3):
                            rows_available_3 = rows_available_2.copy()
                            rows_available_3.remove(row5)
                            rows_available_3.remove(row6)
                            rows_available_3.remove(row7)

                            rows_chosen_3 = rows_chosen_2 + [row5, row6, row7]
                            if not check_row(A, row5, row6, row7, cols_chosen_3, len(cols_chosen_3)/3, rank=desired_rank, check_edges = True):
                                continue

                            for col8, col9, col10 in itertools.combinations(cols_available_3, 3):
                                cols_available_4 = cols_available_3.copy()
                                cols_available_4.remove(col8)
                                cols_available_4.remove(col9)
                                cols_available_4.remove(col10)

                                cols_chosen_4 = cols_chosen_3 + [col8, col9, col10]
                                if not check_column(A, col8, col9, col10, rows_chosen_3, len(rows_chosen_3)/3, rank=desired_rank, check_edges = True):
                                    continue

                                for row8, row9, row10 in itertools.combinations(rows_available_3, 3):
                                    rows_available_4 = rows_available_3.copy()
                                    rows_available_4.remove(row8)
                                    rows_available_4.remove(row9)
                                    rows_available_4.remove(row10)

                                    rows_chosen_4 = rows_chosen_3 + [row8, row9, row10]
                                    if not check_row(A, row8, row9, row10, cols_chosen_4, len(cols_chosen_4)/3, rank=desired_rank, check_edges = True):
                                        continue

                                    for col11, col12, col13 in itertools.combinations(cols_available_4, 3):
                                        cols_available_5 = cols_available_4.copy()
                                        cols_available_5.remove(col11)
                                        cols_available_5.remove(col12)
                                        cols_available_5.remove(col13)

                                        cols_chosen_5 = cols_chosen_4 + [col11, col12, col13]
                                        if not check_column(A, col11, col12, col13, rows_chosen_4, len(rows_chosen_4)/3, rank=desired_rank, check_edges = True):
                                            continue

                                        for row11, row12, row13 in itertools.combinations(rows_available_4, 3):
                                            rows_available_5 = rows_available_4.copy()
                                            rows_available_5.remove(row11)
                                            rows_available_5.remove(row12)
                                            rows_available_5.remove(row13)

                                            rows_chosen_5 = rows_chosen_4 + [row11, row12, row13]
                                            if not check_row(A, row11, row12, row13, cols_chosen_5, len(cols_chosen_5)/3, rank=desired_rank, check_edges = True):
                                                continue

                                            for col14, col15, col16 in itertools.combinations(cols_available_5, 3):
                                                cols_available_6 = cols_available_5.copy()
                                                cols_available_6.remove(col14)
                                                cols_available_6.remove(col15)
                                                cols_available_6.remove(col16)

                                                cols_chosen_6 = cols_chosen_5 + [col14, col15, col16]
                                                if not check_column(A, col14, col15, col16, rows_chosen_5, len(rows_chosen_5)/3, rank=desired_rank, check_edges = True):
                                                    continue

                                                for row14, row15, row16 in itertools.combinations(rows_available_5, 3):
                                                    rows_available_6 = rows_available_5.copy()
                                                    rows_available_6.remove(row14)
                                                    rows_available_6.remove(row15)
                                                    rows_available_6.remove(row16)

                                                    rows_chosen_6 = rows_chosen_5 + [row14, row15, row16]
                                                    if not check_row(A, row14, row15, row16, cols_chosen_6, len(cols_chosen_6)/3, rank=desired_rank, check_edges = True):
                                                        continue

                                                    for col17, col18, col19 in itertools.combinations(cols_available_6, 3):
                                                        cols_available_7 = cols_available_6.copy()
                                                        cols_available_7.remove(col17)
                                                        cols_available_7.remove(col18)
                                                        cols_available_7.remove(col19)

                                                        cols_chosen_7 = cols_chosen_6 + [col17, col18, col19]
                                                        if not check_column(A, col17, col18, col19, rows_chosen_6, len(rows_chosen_6)/3, rank=desired_rank, check_edges = True):
                                                            continue

                                                        for row17, row18, row19 in itertools.combinations(rows_available_6, 3):
                                                            rows_available_7 = rows_available_6.copy()
                                                            rows_available_7.remove(row17)
                                                            rows_available_7.remove(row18)
                                                            rows_available_7.remove(row19)

                                                            rows_chosen_7 = rows_chosen_6 + [row17, row18, row19]
                                                            if not check_row(A, row17, row18, row19, cols_chosen_7, len(cols_chosen_7)/3, rank=desired_rank, check_edges = True):
                                                                continue

                                                            all_cols = [col0, col1] + cols_chosen_7
                                                            all_rows = [row0, row1] + rows_chosen_7
                                                            # perform transformation in place
                                                            return matrix(20, 20, lambda i, j: mat[all_rows[i], all_cols[j]])
    # there is none
    return None

