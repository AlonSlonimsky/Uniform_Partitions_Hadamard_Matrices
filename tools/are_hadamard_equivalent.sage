import subprocess
nauty_location = "../../../nauty2_8_9/dreadnaut"

#------------------------------------------------------
# are_hadamard_equivalent
#
# PURPOSE: Check if two Hadamard matrices are equivalent by using nauty.
# PARAMETERS:
#     Matrix mat1: The first Hadamard matrix.
#     Matrix mat2: The second Hadamard matrix.
# Returns: True if mat1 is Hadamard equivalent to mat2, false otherwise.
#------------------------------------------------------
def are_hadamard_equivalent(mat1: Matrix, mat2: Matrix) -> bool:
    script = write_hadamard_dreadnaut_script(mat1, mat2) # get the neccessary command for nauty
    output = subprocess.run(nauty_location, input = script, shell = True, capture_output = True, text = True) # launch the subprocess and run it
    return output.stdout.find("identical") != -1 # basic check for a substring in the output

#------------------------------------------------------
# write_hadamard_dreadnaut_script
#
# PURPOSE: Get the command to check if two Hadamard matrices are equivalent to 
# PARAMETERS:
#     Matrix mat1: The first Hadamard matrix.
#     Matrix mat2: The second Hadamard matrix.
# Returns: The command, to be fed into nauty.
#------------------------------------------------------
def write_hadamard_dreadnaut_script(mat1: Matrix, mat2: Matrix) -> str:
    def graph_transform(mat):
        n = mat.nrows()

        # map between a vertex and its number
        # must enumerate all vertices to work in nauty
        vertex_map = {
            ("v", i): i for i in range(n)
        }
        vertex_map.update({
            ("v'", i): n+i for i in range(n)
        })
        vertex_map.update({
            ("w", i): 2*n+i for i in range(n)
        })
        vertex_map.update({
            ("w'", i): 3*n+i for i in range(n)
        })

        # set up list of adjacent vertices
        num_vertices = 4*n
        adjacency = [[] for vertex in range(num_vertices)]

        # add loops
        for i in range(n):
            vi = vertex_map[("v", i)]
            vpi = vertex_map[("v'", i)]
            adjacency[vi].append(vi)
            adjacency[vpi].append(vpi)

        # add the rest of the edges
        for i in range(n): # row
            for j in range(n): # col
                if mat[i, j] == 1:
                    adjacency[vertex_map[("v", i)]].append(vertex_map[("w", j)])
                    adjacency[vertex_map[("v'", i)]].append(vertex_map[("w'", j)])
                else:
                    adjacency[vertex_map[("v", i)]].append(vertex_map[("w'", j)])
                    adjacency[vertex_map[("v'", i)]].append(vertex_map[("w", j)])

        # nauty command for this
        lines = [f"n={num_vertices} g"]
        for i in range(num_vertices):
            neighbors = adjacency[i]
            neighbor_str = " ".join(str(j) for j in neighbors) # space between each neighbour
            lines.append(f"{i} : {neighbor_str};")
        return lines

    # flags
    output = ("c -a \n\n")

    # first
    lines1 = graph_transform(mat1)
    output+=("\n".join(lines1)) # newline between each vertex
    output+=("\nx @\n\n")

    # second
    lines2 = graph_transform(mat2)
    output+=("\n".join(lines2)) # newline between each vertex
    output+=("\nx\n\n")

    output+=("##\n")

    return output
