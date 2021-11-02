import numpy as np

def generate_matrix(dim, v1, v2, v3):
    v = np.zeros((3 * dim, dim))
    for i in range(3 * dim):
        if i < dim:
            for j in range(dim):
                v[i, j] = v1
        elif i < 2 * dim:
            for j in range(dim):
                v[i, j] = v2
        else:
            for j in range(dim):
                v[i, j] = v3
    return v

def change_value(matrix, index, i, j, value):
    dim = matrix.shape[1]
    matrix[i + dim * (index - 1) - 1, j - 1] = value
    matrix[j + dim * (index - 1) - 1, i - 1] = value

def write(filename, variable):
    file = open(filename, "w")
    file.write(">>>>>>>>> GB properties <<<<<<<<<\n")
    file.write("---GB energy (sigma_ij), mobility prefactor (mob0_ij), activation energy (Q_ij)---\n")
    shape = variable.shape
    for i in range(shape[0]):
        for j in range(shape[1]):
            if i >= shape[1] and i < 2 * shape[1]:
                file.write(f"{variable[i, j]:.2e} ")
            else:
                file.write(f"{variable[i, j]:.3f} ")
        file.write(f"\n")
    file.close()

if __name__ == "__main__":
    filename = "./poly_anisotropy_mobility.txt"
    v = generate_matrix(15, 0.708, 2.5e-6, 0.23)
    change_value(v, 2, 1, 3, 5.0e-6) # change_value(matrix, index, i, j, value)
    change_value(v, 2, 4, 5, 5.0e-6) # change_value(matrix, index, i, j, value)
    write(filename, v)