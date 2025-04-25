from matrix import Matrix

if __name__ == "__main__":
    
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
