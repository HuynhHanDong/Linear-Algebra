from matrix import *
from vector import *

if __name__ == "__main__":
    # Matrix
    matrixA = Matrix([[1, -2, 2], [-2, 1, -2], [2, -2, 1]])
    matrixA.show_matrix()

    print("Transpose of matrix A:")
    matrixB = matrixA.transpose()
    matrixB.show_matrix()

    print("Inverse of matrix A:")
    matrixC = matrixA.inverse()
    matrixC.show_matrix()

    print("Eigenvalues of matrix A:")
    eigenvalues = matrixA.find_eigenvalues()
    print(eigenvalues)

    print("Eigenvector of matrix A:")
    eigenvector = matrixA.find_eigenvector(eigenvalues[1])
    print(eigenvector)

    print("Diagonalization of matrix A:")
    diagonalize = matrixA.diagonalize()
    for matrix in diagonalize:
        matrix.show_matrix()

    # Vector
    vectorA = Vector([1, -2, 2])
    vectorB = Vector([1, 2, 3])
    dot = vectorA.dot_product(vectorB)
    print("Dot product:", dot)
    print("Cross product:")
    product = vectorA.cross_product(vectorB)
    product.show_vector()