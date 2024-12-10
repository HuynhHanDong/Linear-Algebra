from matrix import Matrix

if __name__ == "__main__":
    '''
    matrixA = Matrix([[2, -3, 3, 2], [3, 0, 1, 4], [-2, 0, 0, 2], [4, 0, -1, 5]])
    matrixA.show_matrix()
    matrixB = Matrix([[1, 2, 3], [3, 2, 1]]).transpose_matrix()
    matrixB.show_matrix()
    matrixC = matrixA.inverse_matrix()
    matrixC.show_matrix() '''
    matrixA = Matrix([[1, 0, 0], [-2, 3, 0], [0, 0, 4]])
    matrixA.show_matrix()
    B = matrixA.find_eigenvector(1)
    print(B)

