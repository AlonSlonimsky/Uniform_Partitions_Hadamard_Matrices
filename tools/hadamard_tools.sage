#------------------------------------------------------
# is_hadamard
#
# PURPOSE: Check if a square +1,-1 matrix is Hadamard.
#          Checks orthogonality via dot products to allow the option to escape early.
# PARAMETERS:
#     Matrix mat: The matrix to be checked. Must be square and only consist of +1s and -1s.
# Returns: True if mat is Hadamard, false otherwise.
#------------------------------------------------------
def is_hadamard(mat: Matrix) -> bool:
    order = mat.nrows()
    for firstRow in range(order):
        for secondRow in range(firstRow+1,order):
            if mat[firstRow].dot_product(mat[secondRow]) != 0:
                return False
    return True

#------------------------------------------------------
# partition_3x3
#
# PURPOSE: Partition a matrix into 3x3 sub-blocks startings from the rightmost column and bottommost row.
# PARAMETERS:
#     Matrix mat: The matrix to be partitioned.
# Returns: The partitioned matrix.
#------------------------------------------------------
def partition_3x3(mat: Matrix) -> Matrix:
    order = mat.nrows()
    subdivisionLines = list(range(order - 3, 0, -3))
    mat.subdivide(subdivisionLines, subdivisionLines)
    return mat

#------------------------------------------------------
# get_3x3_sub_blocks
#
# PURPOSE: Get the 3x3 sub-blocks of the given matrix. 
# PARAMETERS:
#     Matrix mat: The matrix being partitioned.
# Returns: The 3x3 sub-blocks in order left-right, top-bottom. Anything not fitting into a 3x3 block is ignored.
#------------------------------------------------------
def get_3x3_sub_blocks(mat: Matrix) -> list[Matrix]:
    submatrices = []
    for i in range(mat.nrows() - 3, -1, -3):
        for j in range(mat.ncols() - 3, -1, -3):
            sub = mat[i:i+3, j:j+3]
            submatrices.append(sub)
    submatrices.reverse()
    return submatrices

#------------------------------------------------------
# get_matrix_of_ranks
#
# PURPOSE: Get the rank of all the 3x3 sub-blocks of the given matrix as a matrix.
# PARAMETERS:
#     Matrix mat: The matrix being manipulated.
# Returns: A matrix with the ranks of all the 3x3 sub-blocks as elements is corresponding positions. Anything not fitting into a 3x3 block is ignored.
#------------------------------------------------------
def get_matrix_of_3x3_ranks(mat: Matrix) -> Matrix:
    submatrices = get_3x3_sub_blocks(mat)
    n = Integer(len(submatrices)).sqrt()
    return Matrix(ZZ, n, n, [submatrix.rank() for submatrix in submatrices])

#------------------------------------------------------
# get_list_of_ranks
#
# PURPOSE: Get the rank of all the 3x3 sub-blocks of the given matrix as a list.
# PARAMETERS:
#     Matrix mat: The matrix being manipulated.
# Returns: A list with the ranks of all the 3x3 sub-blocks as elements is corresponding positions. Anything not fitting into a 3x3 block is ignored.
#------------------------------------------------------
def get_list_of_3x3_ranks(mat: Matrix) -> list[int]:
    submatrices = get_3x3_sub_blocks(mat)
    n = Integer(len(submatrices)).sqrt()
    return [submatrix.rank() for submatrix in submatrices]


#------------------------------------------------------
# get_matrix_signature
#
# PURPOSE: Get the signature of all the 3x3 sub-blocks of the given matrix as a matrix.
# PARAMETERS:
#     Matrix mat: The matrix being manipulated.
# Returns: A matrix with the signature of all the 3x3 sub-blocks as elements is corresponding positions. Anything not fitting into a 3x3 block is ignored.
#------------------------------------------------------
def get_matrix_signature(mat: Matrix) -> Matrix:
    submatrices = get_3x3_sub_blocks(mat)
    n = Integer(len(submatrices)).sqrt()
    return Matrix(ZZ, n, n, [submatrix.determinant()/4 for submatrix in submatrices])

def get_3x3_quasisignature(mat: Matrix) -> int:
    return  (mat[0,0]*mat[1,1]*mat[2,2]) + (mat[0,1]*mat[1,2]*mat[2,0]) + (mat[0,2]*mat[1,0]*mat[2,1])

def get_matrix_quasisignature(mat: Matrix) -> Matrix:
    submatrices = get_3x3_sub_blocks(mat)
    n = Integer(len(submatrices)).sqrt()
    return Matrix(ZZ, n, n, [get_3x3_quasisignature(submatrix) for submatrix in submatrices])

#------------------------------------------------------
# weave_matrices_basic
#
# PURPOSE: Weave two matrices together where every matrix in the warp is the same, same with the woof.
# PARAMETERS:
#     Matrix mat1: The warp.
#     Matrix mat2: The woof.
# Returns: The result of the basic weaving of the two matrices.
#------------------------------------------------------
def weave_matrices_basic(mat1: Matrix, mat2: Matrix) -> Matrix:
    blocks = []
    for r in range(mat2.nrows()):
        for c in range(mat1.ncols()):
            blocks.append(mat1.submatrix(col = c, ncols = 1) * mat2.submatrix(row = r, nrows = 1))
    return block_matrix(mat2.nrows(), mat1.ncols(), blocks) # block_matrix is a very handy constructor for matrices

#------------------------------------------------------
# kronecker_product
#
# PURPOSE: Return the Kronecker product of two matrices. Actually just use built-in tensor_product function.
# PARAMETERS:
#     Matrix mat1: The first operand.
#     Matrix mat2: The second operand.
# Returns: mat1 ⊗ mat2.
#------------------------------------------------------
def kronecker_product(mat1: Matrix, mat2: Matrix) -> Matrix:
    return mat1.tensor_product(mat2)

#------------------------------------------------------
# normalize_matrix
#
# PURPOSE: Normalize a given matrix by negating rows & columns as needed.
# PARAMETERS:
#     Matrix mat: The matrix to be normalized.
# Returns: The matrix in normalized form.
#------------------------------------------------------
def normalize_matrix(mat: Matrix) -> Matrix:
    for rowIndex in range(mat.nrows()):
        if mat[rowIndex, 0] == -1:
            mat[rowIndex, :] = -mat[rowIndex, :] # negate row
    
    
    for colIndex in range(mat.ncols()):
        if mat[0, colIndex] == -1:
            mat[:, colIndex] = -mat[:, colIndex] # negate column
        
    #technically an extra redundant check here on the (0,0)th position, but not worth caring about

    return mat

#------------------------------------------------------
# compress_cols
#
# PURPOSE: Compress the columns of a matrix by their type:
# (1,1,1) = a = (-1,-1,-1)
# (-1,1,1) = b = (1,-1,-1)
# (1,-1,1) = c = (-1,1,-1)
# (1,1,-1) = d = (-1,-1,1)
# NOTE: this is a 2d array of chars, not a sagemath matrix
# PARAMETERS:
#     Matrix mat: The matrix to be compressed.
# Returns: The compressed matrix
#------------------------------------------------------
def compress_cols(mat: Matrix):
    #will compress by actually just rebuilding the matrix from scratch
    a = matrix(3,1,[-1,-1,-1])
    b = matrix(3,1,[-1,1,1])
    c = matrix(3,1,[1,-1,1])
    d = matrix(3,1,[1,1,-1])

    shift = mat.ncols() % 3

    result = []
    for r in range(floor(mat.nrows()/3)):
        for co in range(mat.ncols()):
            col = mat.submatrix(row = r*3+shift, col = co, nrows = 3, ncols = 1)
            if col == a:
                result.append('a')
            elif col == b:
                result.append('b')
            elif col == c:
                result.append('c')
            elif col == d:
                result.append('d')
            elif col == -a:
                result.append('a')
            elif col == -b:
                result.append('b')
            elif col == -c:
                result.append('c')
            elif col == -d:
                result.append('d')
    
    n = mat.ncols()
    num_rows = n//3

    return [result[i * n:(i + 1) * n] for i in range(num_rows)]


#------------------------------------------------------
# load_partitioned
#
# PURPOSE: Load a partitioned matrix from a string
# PARAMETERS:
#     str string: The string containing the partitioned matrix.
# Returns: The matrix in normalized form.
#------------------------------------------------------
def load_partitioned(string: str, order = 24) -> Matrix:
    element_list = []
    string = string.replace("|", " ").replace("[", " ").replace("]", " ") # make the entire string whitespace apart from the numbers
    string = string.split()
    for char in string:
        if char == "1":
            element_list.append(1)
        elif char == "-1":
            element_list.append(-1)
    return Matrix(order,order,element_list)
