class Matrix:
    def __init__(self, matrix):
        self.matrix = matrix
        self.rows = len(matrix)
        self.cols = len(matrix[0])

    def dimension(self):
        return self.rows, self.cols

    def show_matrix(self):
        for i in range(self.rows):
            for j in range(self.cols):
                print(self.matrix[i][j], end=' ')
            print()
        print()

    def transpose(self):
        transpose = [[self.matrix[i][j] for i in range(self.rows)] for j in range(self.cols)]
        return Matrix(transpose)

    def add(self, other):
        if self.dimension() != other.dimension():
            print("Matrices must have the same dimensions")
            return
        result = [[self.matrix[i][j] + other.matrix[i][j] for j in range(self.cols)] for i in range(self.rows)]
        return Matrix(result)

    def subtract(self, other):
        if self.dimension() != other.dimension():
            print("Matrices must have the same dimensions")
            return
        result = [[self.matrix[i][j] - other.matrix[i][j] for j in range(self.cols)] for i in range(self.rows)]
        return Matrix(result)

    def scalar_multiple(self, k):
        result = [[self.matrix[i][j] * k for j in range(self.cols)] for i in range(self.rows)]
        return Matrix(result)

    def multiply(self, other):
        if self.cols != other.rows:
            print("Columns in first matrix must equal to number of rows in second matrix")
            return
        result = [[sum(self.matrix[i][k] * other.matrix[k][j] for k in range(self.cols)) for j in range(other.cols)] for i in range(self.rows)]
        return Matrix(result)

    @staticmethod
    def identity(n):
        identity = [[1 if i == j else 0 for i in range(n)] for j in range(n)]
        return Matrix(identity)

    def gaussian_elimination(self):
        mat = self.matrix
        n = min(self.dimension())
        for i in range(n):
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
        return mat

    def rank(self):
        rank = 0
        mat = self.gaussian_elimination()
        for i in range(len(mat)):
            if any(mat[i][j] != 0 for j in range(len(mat[0]))):
                rank += 1
        return rank

    def is_square(self):
        if self.rows == self.cols:
            return True
        else:
            return False

    def determinant(self):
        if not self.is_square():
            print("Matrix must be square")
            return
        mat = self.matrix
        n = len(mat)
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

    def inverse(self):
        if not self.is_square():
            print("Matrix must be square")
            return
        mat = self.matrix
        n = len(mat)
        identity = [[1 if i == j else 0 for i in range(n)] for j in range(n)]

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
                    print("Det = 0 -> Not invertible!")
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
        inverse = [[round(element, 4) for element in row[n:]] for row in mat]
        return Matrix(inverse)

    ''' gram_schmidt chưa tổng quát'''
    def gram_schmidt(self):
        n = self.rows
        A = [[0] * n for _ in range(n)]
        B = [[0] * n for _ in range(n)]
        for j in range(n):
            v = Vector([self.matrix[i][j] for i in range(n)])
            for i in range(j):
                a = Vector([A[k][i] for k in range(n)])
                B[i][j] = a.dot_product(v)
                v = v.vector_subtract(a.scalar_multiply(B[i][j]))

            B[j][j] = v.norm_vector(2)
            rounded = [round(x / B[j][j]) for x in v]

            for i in range(n):
                A[i][j] = rounded[i]

        return A, B

    def qr_algorithm(self, num_iterations):
        n = self.rows
        A = self.matrix

        for _ in range(num_iterations):
            # Phân rã ma trận A thành QR bằng Gram-Schmidt
            Q, R = self.gram_schmidt()
            A = [[sum(R[i][k] * Q[k][j] for k in range(n)) for j in range(n)] for i in range(n)]

        return A

    def find_eigenvalue(self, num_iterations=100):
        n = self.rows
        A = self.qr_algorithm(num_iterations)
        # Các giá trị trên đường chéo của A là eigenvalues
        eigenvalues = [A[i][i] for i in range(n)]
        return eigenvalues

    def find_eigenvector(self, eigenvalue):
        # Tìm vector riêng tương ứng với một giá trị riêng
        n = self.rows
        A = self.matrix

        # Tạo ma trận A - λI
        B = [[A[i][j] - (eigenvalue if i == j else 0) for j in range(n)] for i in range(n)]

        # Giả sử giá trị riêng đã là đúng, ta tìm vector riêng cho ma trận A - λI
        # Bằng cách giải hệ phương trình Bx = 0
        # Ở đây ta dùng phương pháp khử Gauss để tìm x
        B = Matrix(B).gaussian_elimination()

        x = [0 for _ in range(n)]
        x[-1] = 1  # Ta giả sử giá trị cuối cùng của vector riêng là 1
        for i in range(n-2, -1, -1):
            x[i] = -sum(B[i][j] * x[j] for j in range(i+1, n)) / B[i][i]
        eigenvector = x

        return eigenvector

    def diagonalize(self, num_iterations=100):
        # Tính giá trị riêng
        eigenvalues = self.find_eigenvalue(num_iterations)

        # Tìm vector riêng tương ứng cho mỗi giá trị riêng
        eigenvectors = []
        for eigenvalue in eigenvalues:
            eigenvector = self.find_eigenvector(eigenvalue)
            eigenvectors.append(eigenvector)

        # Ma trận chéo hóa P, chứa các vector riêng
        P = [[eigenvectors[j][i] for j in range(len(eigenvalues))] for i in range(len(eigenvalues))]

        # Tạo bản sao của P
        P_inv = Matrix(P).inverse()

        # Ma trận chéo chứa các giá trị riêng
        Lambda = [[eigenvalues[i] if i == j else 0 for j in range(len(eigenvalues))] for i in range(len(eigenvalues))]

        return Matrix(P), Matrix(Lambda), P_inv

    def is_symmetric(self):
        if self == self.transpose():
            return True
        else:
            return False

    def is_positive_definite(self):
        if self.is_square() and self.is_symmetric():
            return True
        return False

    def norm_matrix(self):
        # Norm Frobenius
        return (sum(cell ** 2 for row in self.matrix for cell in row)) ** 0.5

class Vector:
    def __init__(self, vector):
        self.vector = vector
        self.size = len(vector)

    def dot_product(self, other):
        dotProduct = 0
        for i in range(self.size):
            dotProduct += self.vector[i] * other.vector[i]
        return dotProduct

    def scalar_multiply(self, k):
        return [k * x for x in self.vector]

    def vector_subtract(self, other):
        ans = [0 for _ in range(self.size)]
        for i in range(self.size):
            ans[i] += self.vector[i] - other.vector[i]
        return ans

    def norm_vector(self, p=2):
        return (sum(abs(x)**p for x in self.vector))**(1/p)
