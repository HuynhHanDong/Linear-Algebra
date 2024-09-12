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
                    print("Det = 0. Not invertible!")
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

    '''
    def gram_schmidt(self, A):
        """Phân rã ma trận A thành QR bằng Gram-Schmidt"""
        n = len(A)
        Q = [[0] * n for _ in range(n)]
        R = [[0] * n for _ in range(n)]

        for j in range(n):
            v = [A[i][j] for i in range(n)]

            for i in range(j):
                q = [Q[k][i] for k in range(n)]
                R[i][j] = dot_product(q, v)
                v = self.subtract(v, self.scalar_multiple(q, R[i][j]))

            R[j][j] = norm_vector(v)
            q = [x / R[j][j] for x in v]

            for i in range(n):
                Q[i][j] = q[i]

        return Q, R

    def qr_algorithm(self, A, num_iterations):
        """Sử dụng thuật toán QR để tìm giá trị riêng"""
        n = len(A)

        for _ in range(num_iterations):
            Q, R = self.gram_schmidt(A)
            A = [[sum(R[i][k] * Q[k][j] for k in range(n)) for j in range(n)] for i in range(n)]

        # Các giá trị trên đường chéo của A là eigenvalues
        eigenvalues = [A[i][i] for i in range(n)]
        return eigenvalues, A

    def find_eigenvector(self, A, eigenvalue):
        """Tìm vector riêng tương ứng với một giá trị riêng"""
        n = len(A)
        I = [[1 if i == j else 0 for j in range(n)] for i in range(n)]

        # Tạo ma trận A - λI
        B = [[A[i][j] - (eigenvalue if i == j else 0) for j in range(n)] for i in range(n)]

        # Giả sử giá trị riêng đã là đúng, ta tìm vector riêng cho ma trận A - λI
        # Bằng cách giải hệ phương trình Bx = 0
        # Ở đây ta dùng phương pháp khử Gauss để tìm x

    def diagonalize(self, A, num_iterations=100):
        """Tính giá trị riêng và chéo hóa ma trận"""
        # Tính giá trị riêng bằng QR algorithm
        eigenvalues, _ = self.qr_algorithm(A, num_iterations)

        # Tìm vector riêng tương ứng cho mỗi giá trị riêng
        eigenvectors = []
        for eigenvalue in eigenvalues:
            eigenvector = self.find_eigenvector(A, eigenvalue)
            eigenvectors.append(eigenvector)

        # Ma trận chéo hóa P, chứa các vector riêng
        P = [[eigenvectors[j][i] for j in range(len(eigenvalues))] for i in range(len(eigenvalues))]

        P_inv = self.inverse_matrix([row[:] for row in P])  # Tạo bản sao của P

        # Ma trận chéo chứa các giá trị riêng
        Lambda = [[eigenvalues[i] if i == j else 0 for j in range(len(eigenvalues))] for i in range(len(eigenvalues))]

        return P, Lambda, P_inv
    '''

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

def dot_product(vectorA, vectorB):
    n = len(vectorA)
    dotProduct = 0
    for i in range(n):
        dotProduct += vectorA[i] * vectorB[i]
    return dotProduct

def scalar_multiply(vector, k):
    return [k * x for x in vector]

def vector_subtract(vectorA, vectorB):
    n = len(vectorA)
    ans = 0
    for i in range(n):
        ans += vectorA[i] - vectorB[i]
    return ans

def norm_vector(vector, p):
    return (sum(abs(x)**p for x in vector))**(1/p)
