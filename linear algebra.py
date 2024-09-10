def create_matrix():
    row = int(input("Enter number of row: "))
    col = int(input("Enter number of col: "))
    matrix = [[0 for i in range(col)] for j in range(row)]
    for i in range(row):
        for j in range(col):
            matrix[i][j] = int(input("Enter value: "))
    return matrix

def matrix_dimension(matrix):
    row = len(matrix)
    col = len(matrix[0])
    return row, col

def show_matrix(matrix):
    row, col = matrix_dimension(matrix)
    for i in range(col):
        for j in range(row):
            print(matrix[i][j], end=' ')
        print()
    print()

def transpose_matrix(matrix):
    row, col = matrix_dimension(matrix)
    transpose = [row[:] for row in matrix]
    for i in range(row):
        for j in range(col):
            transpose[i][j] = transpose[j][i]
    return transpose

def add(matrixA, matrixB):
    rowA, colA = matrix_dimension(matrixA)
    rowB, colB = matrix_dimension(matrixB)
    if rowA != rowB or colA != colB:
        print("2 matrices don't the same dimension")
        return
    new = [[0 for i in range(colA)] for j in range(rowA)]
    for i in range(rowA):
        for j in range(colA):
            new[i][j] = matrixA[i][j] + matrixB
    return new

def subtract(matrixA, matrixB):
    rowA, colA = matrix_dimension(matrixA)
    rowB, colB = matrix_dimension(matrixB)
    if rowA != rowB or colA != colB:
        print("2 matrices don't the same dimension")
        return
    new = [[0 for i in range(colA)] for j in range(rowA)]
    for i in range(rowA):
        for j in range(colA):
            new[i][j] = matrixA[i][j] - matrixB
    return new

def scalar_multiple(matrix, k):
    row, col = matrix_dimension(matrix)
    new = [row[:] for row in matrix]
    for i in range(row):
        for j in range(col):
            new[i][j] *= k
    return new

def product_matrix(matrixA, matrixB):
    rowA, colA = matrix_dimension(matrixA)
    rowB, colB = matrix_dimension(matrixB)
    if colA != rowB:
        print("Cannot multiply")
        return
    new = [[0 for i in range(colA)] for j in range(rowA)]
    c = 0
    for i in range(rowA):
        for j in range(colB):
            for k in range(colA):
                k += matrixA[i][k] * matrixB[k][j]
            new[i][j] = c
            c = 0
    return new

def identity_matrix(n):
    identity = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        identity[i][i] = 1
    return identity

def determinant(matrix):
    n = len(matrix)
    mat = [row[:] for row in matrix]
    det = 1

    # Gaussian elimination
    for i in range(n):
        if mat[i][i] == 0:
            for j in range(i + 1, n):
                if mat[j][i] != 0:
                    mat[i], mat[j] = mat[j], mat[i]
                    det *= -1
                    break
            else:
                return 0

        det *= mat[i][i]

        for j in range(i + 1, n):
            ratio = mat[j][i] / mat[i][i]
            for k in range(i, n):
                mat[j][k] -= ratio * mat[i][k]

    return round(det)

def rank(m, n, matrix):
    mat = [r[:] for r in matrix]
    rank = 0

    # Gaussian elimination
    for i in range(min(n, m)):
        if mat[i][i] == 0:
            for j in range(i + 1, n):
                if mat[j][i] != 0:
                    mat[i], mat[j] = mat[j], mat[i]
                    break
        if mat[i][i] == 0:
            continue

        for j in range(i + 1, n):
            ratio = mat[j][i] / mat[i][i]
            for k in range(i, n):
                mat[j][k] -= ratio * mat[i][k]

    for i in range(m):
        if any(mat[i][j] != 0 for j in range(n)):
            rank += 1

    return rank

def inverse_matrix(matrix):
    n = len(matrix)
    mat = [row[:] for row in matrix]

    # Tạo ma trận đơn vị cùng kích thước
    identity = [[1 if i == j else 0 for j in range(n)] for i in range(n)]

    # Kết hợp ma trận ban đầu và ma trận đơn vị
    for i in range(n):
        mat[i] += identity[i]

    # Biến đổi Gauss-Jordan
    for i in range(n):
        if mat[i][i] == 0:
            for j in range(i + 1, n):
                if mat[j][i] != 0:
                    mat[i], mat[j] = mat[j], mat[i]
                    break
            else:
                print("Not invertible!")
                return

        # Chia dòng hiện tại cho phần tử chính để phần tử đó bằng 1
        pivot = mat[i][i]
        for j in range(2 * n):
            mat[i][j] /= pivot

        # Khử tất cả các phần tử khác trên và dưới phần tử chính
        for j in range(n):
            if j != i:
                ratio = mat[j][i]
                for k in range(2 * n):
                    mat[j][k] -= ratio * mat[i][k]

    # Tách ma trận nghịch đảo ra khỏi ma trận mở rộng
    inverse = [[round(element, 0) for element in row[n:]] for row in mat]
    return inverse

def eigenvalue(m, n, a):
    pass

def diagonalization(m, n, a):
    pass

def positive_definite(m, n, a):
    pass

def norm_matrix(matrix):
    # Norm Frobenius
    return (sum(cell ** 2 for row in matrix for cell in row))**0.5

def norm_vector(vector, p):
    return (sum(abs(x)**p for x in vector))**(1/p)

def main():
    row = 2
    col = 3
    matrixA = [[2, -3, 3, 2], [3, 0, 1, 4], [-2, 0, 0, 2], [4, 0, -1, 5]]
    row2 = 3
    col2 = 4
    matrixB = [[3, 1], [5, 2]]
    show_matrix(matrixA)

'''
show_matrix(transpose_matrix(matrix))
show_matrix(product_matrix(matrixA, matrixB))
print(subtract(matrixA, matrixB))
print(scalar_multiple(matrixA, 5)
show_matrix(identity_matrix(4))
print("Det = ",determinant(matrix))
print("Rank = ", rank(matrix))
show_matrix(inverse_matrix(matrix))
'''

if __name__ == "__main__":
    main()
