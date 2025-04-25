from vector import *

class Matrix:
    """
    A class representing a matrix with various linear algebra operations.
    
    Attributes:
        matrix (list): A 2D list representing the matrix elements.
        rows (int): Number of rows in the matrix.
        cols (int): Number of columns in the matrix.
    """
    def __init__(self, matrix: list):
        """
        Initialize a matrix object from a 2D list.
        
        Parameters:
            matrix (list): A 2D list representing the matrix elements.
        """
        self.matrix = matrix
        self.rows = len(matrix)
        self.cols = len(matrix[0])

    def dimension(self):
        """
        Get the dimensions of the matrix.
        
        Returns:
            tuple: A tuple containing the number of rows and columns (rows, cols).
        """
        return self.rows, self.cols

    def show_matrix(self):
        """
        Display the matrix elements.
        Prints each row of the matrix on a separate line.
        """
        for i in range(self.rows):
            for j in range(self.cols):
                print(self.matrix[i][j], end=' ')
            print()
        print()

    def __eq__(self, other):
        """
        Check if two matrices are equal.
        
        Parameters:
            other (Matrix): The matrix to compare with.
            
        Returns:
            bool: True if matrices are equal, False otherwise.
        """
        if self.rows != other.rows or self.cols != other.cols:
            return False
        
        for i in range(self.rows):
            for j in range(self.cols):
                if abs(self.matrix[i][j] - other.matrix[i][j]) > 1e-10:
                    return False
        return True

    def transpose(self):
        """
        Calculate the transpose of the matrix.
        
        Returns:
            Matrix: The transposed matrix.
        """
        transpose = [[self.matrix[i][j] for i in range(self.rows)] for j in range(self.cols)]
        return Matrix(transpose)

    def add(self, other):
        """
        Add another matrix to this matrix.
        
        Parameters:
            other (Matrix): The matrix to add.
            
        Returns:
            Matrix: The resulting matrix after addition.
            
        Raises:
            ValueError: If the matrices have different dimensions.
        """
        if self.dimension() != other.dimension():
            raise ValueError("Matrices must have the same dimensions")
        
        result = [[self.matrix[i][j] + other.matrix[i][j] for j in range(self.cols)] for i in range(self.rows)]
        return Matrix(result)

    def subtract(self, other):
        """
        Subtract another matrix from this matrix.
        
        Parameters:
            other (Matrix): The matrix to subtract from this matrix.
            
        Returns:
            Matrix: The resulting matrix after subtraction.
            
        Raises:
            ValueError: If the matrices have different dimensions.
        """
        if self.dimension() != other.dimension():
            raise ValueError("Matrices must have the same dimensions")
        result = [[self.matrix[i][j] - other.matrix[i][j] for j in range(self.cols)] for i in range(self.rows)]
        return Matrix(result)

    def scalar_multiple(self, k: int|float):
        """
        Multiply the matrix by a scalar value.
        
        Parameters:
            k (int | float): The scalar value to multiply the matrix by.
            
        Returns:
            Matrix: A new Matrix object representing the scaled matrix.
        """
        result = [[self.matrix[i][j] * k for j in range(self.cols)] for i in range(self.rows)]
        return Matrix(result)

    def multiply(self, other):
        """
        Multiply this matrix by another matrix.
        
        Parameters:
            other (Matrix): The matrix to multiply with this matrix.
            
        Returns:
            Matrix: A new Matrix object representing the product.
            
        Raises:
            ValueError: If the number of columns in this matrix doesn't equal the number of rows in the other matrix.
        """
        if self.cols != other.rows:
            raise ValueError("Columns in first matrix must equal to number of rows in second matrix")
        
        result = [[sum(self.matrix[i][k] * other.matrix[k][j] for k in range(self.cols)) for j in range(other.cols)] for i in range(self.rows)]
        return Matrix(result)

    @staticmethod
    def identity(n: int):
        """
        Create an identity matrix of size n×n.
        
        Parameters:
            n (int): The size of the identity matrix.
            
        Returns:
            Matrix: A new n×n Matrix object with 1s on the diagonal and 0s elsewhere.
        """
        identity = [[1 if i == j else 0 for i in range(n)] for j in range(n)]
        return Matrix(identity)

    def gaussian_elimination(self, augmented=False):
        """
        Perform Gaussian elimination on the matrix.
        
        Parameters:
            augmented (bool, default=False): If True, treats the matrix as an augmented matrix (system of equations).
            
        Returns:
            Matrix: The matrix in row echelon form.
        """
        # Create a deep copy of the matrix to avoid modifying the original
        mat = [row[:] for row in self.matrix]
        rows, cols = self.rows, self.cols
        
        # Keep track of row swaps for determinant calculation
        row_swaps = 0
        
        # Forward elimination
        pivot_row = 0
        for col in range(cols if not augmented else cols - 1):
            # Find the pivot (maximum element in the current column)
            max_row = pivot_row
            max_val = abs(mat[pivot_row][col]) if pivot_row < rows else 0
            
            for i in range(pivot_row + 1, rows):
                if abs(mat[i][col]) > max_val:
                    max_row = i
                    max_val = abs(mat[i][col])
            
            # If the maximum value is zero, the column is already eliminated
            if abs(max_val) < 1e-10:
                continue
            
            # Swap rows if needed
            if max_row != pivot_row:
                mat[pivot_row], mat[max_row] = mat[max_row], mat[pivot_row]
                row_swaps += 1
            
            # Eliminate all rows below
            for i in range(pivot_row + 1, rows):
                if abs(mat[i][col]) < 1e-10:
                    continue
                    
                factor = mat[i][col] / mat[pivot_row][col]
                mat[i][col] = 0  # Set explicitly to zero (avoid floating point issues)
                
                for j in range(col + 1, cols):
                    mat[i][j] -= factor * mat[pivot_row][j]
            
            pivot_row += 1
            if pivot_row >= rows:
                break
    
        return Matrix(mat)

    def rank(self):
        """
        Calculate the rank of the matrix.
        
        Returns:
            int: The rank of the matrix (number of linearly independent rows/columns).
        """
        rank = 0
        mat = self.gaussian_elimination()
        for i in range(len(mat)):
            if any(mat[i][j] != 0 for j in range(len(mat[0]))):
                rank += 1
        return rank

    def is_square(self):
        """
        Check if the matrix is square (same number of rows and columns).
        
        Returns:
            bool: True if the matrix is square, False otherwise.
        """
        if self.rows == self.cols:
            return True
        else:
            return False

    def determinant(self) -> float:
        """
        Calculate the determinant of the matrix.
        
        Returns:
            float: The determinant of the matrix.
            
        Raises:
            ValueError: If the matrix is not square.
        """
        if not self.is_square():
            raise ValueError("Determinant is defined only for square matrices")
        
        # Create a copy of the matrix
        mat = [row[:] for row in self.matrix]
        n = self.rows
        
        # Keep track of row swaps
        row_swaps = 0
        
        # Perform Gaussian elimination directly in this method
        for i in range(n):
            # Find maximum element in this column
            max_row = i
            for j in range(i + 1, n):
                if abs(mat[j][i]) > abs(mat[max_row][i]):
                    max_row = j
            
            # Swap rows if needed
            if max_row != i:
                mat[i], mat[max_row] = mat[max_row], mat[i]
                row_swaps += 1
            
            # If the pivot is zero, determinant is zero
            if abs(mat[i][i]) < 1e-10:
                return 0
            
            # Eliminate rows below
            for j in range(i + 1, n):
                if abs(mat[i][i]) < 1e-10:
                    continue  # Skip division by near-zero
                factor = mat[j][i] / mat[i][i]
                mat[j][i] = 0  # Set explicitly to avoid floating-point errors
                for k in range(i + 1, n):
                    mat[j][k] -= factor * mat[i][k]
        
        # Calculate determinant as product of diagonal elements
        det = 1.0
        for i in range(n):
            det *= mat[i][i]
        
        # Account for row swaps
        det *= (-1) ** row_swaps
        
        # Handle floating point errors for integer determinants
        if abs(det - round(det)) < 1e-10:
            return round(det)
        return det

    def inverse(self):
        """
        Calculate the inverse of the matrix.
        
        Returns:
            Matrix: The inverse of the matrix
            
        Raises:
            ValueError: If the matrix is not square or not invertible.
        """
        if not self.is_square():
            raise ValueError("Matrix must be square to have an inverse")
            
        n = self.rows
        
        # Create an augmented matrix [A|I]
        augmented = []
        for i in range(n):
            row = self.matrix[i].copy()
            # Append identity matrix
            for j in range(n):
                row.append(1.0 if i == j else 0.0)
            augmented.append(row)
        
        augmented_matrix = Matrix(augmented)
        
        # Forward elimination
        pivot_row = 0
        for col in range(n):
            # Find pivot
            max_row = pivot_row
            max_val = abs(augmented_matrix.matrix[pivot_row][col])
            
            for i in range(pivot_row + 1, n):
                if abs(augmented_matrix.matrix[i][col]) > max_val:
                    max_row = i
                    max_val = abs(augmented_matrix.matrix[i][col])
            
            if abs(max_val) < 1e-10:
                raise ValueError("Matrix is singular and cannot be inverted")
            
            # Swap rows if needed
            if max_row != pivot_row:
                augmented_matrix.matrix[pivot_row], augmented_matrix.matrix[max_row] = \
                    augmented_matrix.matrix[max_row], augmented_matrix.matrix[pivot_row]
            
            # Scale pivot row to make pivot = 1
            pivot = augmented_matrix.matrix[pivot_row][col]
            for j in range(col, 2*n):
                augmented_matrix.matrix[pivot_row][j] /= pivot
            
            # Eliminate in all other rows
            for i in range(n):
                if i != pivot_row:
                    factor = augmented_matrix.matrix[i][col]
                    for j in range(col, 2*n):
                        augmented_matrix.matrix[i][j] -= factor * augmented_matrix.matrix[pivot_row][j]
            
            pivot_row += 1
        
        # Extract the inverse from the right side of the augmented matrix
        inverse = [[augmented_matrix.matrix[i][j+n] for j in range(n)] for i in range(n)]
        
        # Clean up near-zero elements
        for i in range(n):
            for j in range(n):
                if abs(inverse[i][j]) < 1e-10:
                    inverse[i][j] = 0
                elif abs(inverse[i][j] - round(inverse[i][j])) < 1e-10:
                    inverse[i][j] = round(inverse[i][j])
        
        return Matrix(inverse)

    def gram_schmidt(self):
        """
        Perform the Gram-Schmidt orthogonalization process.
        
        Returns:
            tuple: A tuple (Q, R) where Q is an orthogonal matrix and R is upper triangular.
            
        Raises:
            ValueError: If the matrix is not square or its columns are linearly dependent.
        """
        if not self.is_square():
            raise ValueError("Gram-Schmidt process requires a square matrix")
        
        n = self.rows

        # Extract columns as vectors
        columns = [[self.matrix[i][j] for i in range(n)] for j in range(n)]

        # Initialize Q and R matrices
        Q = [[0] * n for _ in range(n)]
        R = [[0] * n for _ in range(n)]

        # Orthonormal basis vectors
        orthonormal_basis = []

        # Apply Gram-Schmidt process
        for j in range(n):
            # Start with the original column vector
            v = Vector(columns[j])

            # Subtract projections of previous orthogonal vectors
            for i in range(j):
                # Calculate projection coefficient (dot product)
                R[i][j] = orthonormal_basis[i].dot_product(v)

                # Subtract projection from current vector
                projection = orthonormal_basis[i].scalar_multiply(R[i][j])
                v = v.subtract(projection)

            # Calculate the norm of the resulting vector
            norm = v.norm(p=2)
            
            if abs(norm) < 1e-10:  # Handle linear dependence
                raise ValueError("Matrix columns are linearly dependent")
            
            # Store the diagonal element of R
            R[j][j] = norm
            
            # Create a normalized vector and store it
            q_j = v.scalar_multiply(1/norm)
            orthonormal_basis.append(q_j)
            
            # Store in Q matrix
            for i in range(n):
                Q[i][j] = q_j.vector[i]

        return Matrix(Q), Matrix(R)

    def qr_algorithm(self, num_iterations, tol=1e-10):
        """
        Apply the QR algorithm to find eigenvalues.

        Parameters:
            num_iterations (int): Maximum number of iterations to perform.
            tol (float, default=1e-10): Tolerance for convergence.
        
        Returns:
            Matrix: A matrix that approximates a diagonal matrix containing eigenvalues.
            
        Raises:
            ValueError: If the matrix is not square.
        """
        if not self.is_square():
            raise ValueError("QR algorithm requires a square matrix")
        
        n = self.rows
        A = [row[:] for row in self.matrix]

        for _ in range(num_iterations):
            # Create a Matrix object from current A
            A_matrix = Matrix(A)

            # Compute QR decomposition
            Q, R = A_matrix.gram_schmidt()
            
            # Compute new A = RQ
            A_new = [[0 for _ in range(n)] for _ in range(n)]
            for i in range(n):
                for j in range(n):
                    A_new[i][j] = sum(R.matrix[i][k] * Q.matrix[k][j] for k in range(n))
            
            # Check for convergence (off-diagonal elements approaching zero)
            if all(abs(A_new[i][j]) < tol for i in range(n) for j in range(n) if i != j):
                A = A_new
                break
                
            A = A_new
        
        return Matrix(A)

    def find_eigenvalues(self, num_iterations=100, tol=1e-10):
        """
        Find eigenvalues of the matrix using the QR algorithm.
        
        Parameters:
            num_iterations (int, default=100): Maximum number of iterations for the QR algorithm.
            tol (float, default=1e-10); Tolerance for convergence.
            
        Returns:
            list: List of eigenvalues.
            
        Raises:
        ValueError: If the matrix is not square.
        """
        if not self.is_square():
            raise ValueError("Eigenvalue calculation requires a square matrix")
            
        # Apply QR algorithm to find Schur form (upper triangular)
        result = self.qr_algorithm(num_iterations, tol)
        
        # Extract eigenvalues from the diagonal
        eigenvalues = [result.matrix[i][i] for i in range(self.rows)]
        
        return eigenvalues

    def find_eigenvector(self, eigenvalue: int|float, tol=1e-10):
        """
        Find an eigenvector corresponding to the given eigenvalue
        using inverse iteration method.
        
        Parameters:
            eigenvalue (int | float): The eigenvalue for which to find the corresponding eigenvector
            tol (float, default=1e-10): Tolerance for numerical computations.
            
        Returns:
            list: The eigenvector corresponding to the given eigenvalue.
            
        Raises:
            ValueError: If the matrix is not square.
        """
        if not self.is_square():
            raise ValueError("Eigenvector calculation requires a square matrix")
            
        n = self.rows
        
        # Create matrix A - λI
        B = [[self.matrix[i][j] - (eigenvalue if i == j else 0) for j in range(n)] for i in range(n)]
        
        # Start with a random vector (ones vector)
        x = [1] * n  
        
        # Normalize
        norm = (sum(val**2 for val in x))**0.5
        x = [val/norm for val in x]
        
        # Inverse iteration
        for _ in range(20):  # Usually converges quickly
            # Solve (A - λI)y = x approximately using Gaussian elimination
            # For simplicity, we'll use a slightly shifted eigenvalue to avoid singularity
            shifted_B = [[B[i][j] + (tol if i == j else 0) for j in range(n)] for i in range(n)]
            
            # Use Gaussian elimination to solve the system
            # This is a simplified approach; a production system would use more robust methods
            augmented = [row[:] + [x[i]] for i, row in enumerate(shifted_B)]
            augmented_matrix = Matrix(augmented)
            solved = augmented_matrix.gaussian_elimination()
            
            # Extract solution y
            y = [0] * n  

            # Perform back-substitution in reverse order
            for i in range(n-1, -1, -1):  
                if abs(solved.matrix[i][i]) < tol:  # Near-zero pivot
                    y[i] = 0
                else:
                    # Calculate the sum of already known values
                    sum_val = 0
                    for j in range(i+1, n):
                        sum_val += solved.matrix[i][j] * y[j]
                    
                    y[i] = (solved.matrix[i][-1] - sum_val) / solved.matrix[i][i]
            
            # Normalize y
            norm = (sum(val**2 for val in y))**0.5
            if norm < tol:  # If norm is too small, retry with different vector
                x = [1 if i == 0 else 0 for i in range(n)]
                continue
                
            x = [val/norm for val in y]
        
        return x

    def diagonalize(self, num_iterations=100, tol=1e-10):
        """
        Attempt to diagonalize the matrix if possible.
        
        Parameters:
            num_iterations (int, default=100): Maximum number of iterations for eigenvalue calculations.
            tol (float, default=1e-10): Tolerance for numerical computations.
            
        Returns:
            tuple: A tuple (P, D, P_inv) where P is the eigenvector matrix,
                  D is the diagonal matrix of eigenvalues, and P_inv is the inverse of P.
            
        Raises:
            ValueError: If the matrix is not square or not diagonalizable.
        """
        if not self.is_square():
            raise ValueError("Only square matrices can be diagonalized")
            
        n = self.rows
        
        # Get eigenvalues
        eigenvalues = self.find_eigenvalues(num_iterations, tol)
        
        # Get eigenvectors for each eigenvalue
        eigenvectors = []
        for eigenvalue in eigenvalues:
            eigenvector = self.find_eigenvector(eigenvalue, tol)
            
            # Normalize eigenvector
            norm = sum(x**2 for x in eigenvector) ** 0.5
            eigenvector = [x/norm for x in eigenvector]
            
            eigenvectors.append(eigenvector)
        
        # Create matrix P (eigenvectors as columns)
        P = [[eigenvectors[j][i] for j in range(n)] for i in range(n)]
        P_matrix = Matrix(P)
        
        # Create diagonal matrix D (eigenvalues on diagonal)
        D = [[eigenvalues[i] if i == j else 0 for j in range(n)] for i in range(n)]
        D_matrix = Matrix(D)
        
        # Calculate inverse of P
        try:
            P_inv = P_matrix.inverse()
        except ValueError:
            raise ValueError("Matrix is not diagonalizable - eigenvectors are linearly dependent")
        
        # Verify diagonalization: P⁻¹AP = D
        # This check helps catch numerical issues
        A_reconstructed = P_inv.multiply(self).multiply(P_matrix)
        if not all(abs(A_reconstructed.matrix[i][j] - D_matrix.matrix[i][j]) < tol 
                for i in range(n) for j in range(n)):
            raise ValueError("Diagonalization failed due to numerical issues")
        
        return P_matrix, D_matrix, P_inv

    def is_symmetric(self):
        """
        Check if the matrix is symmetric (equal to its transpose).
        
        Returns:
            bool: True if the matrix is symmetric, False otherwise.
        """
        if self == self.transpose():
            return True
        else:
            return False
    
    def is_positive_definite(self):
        """
        Check if the matrix is positive definite.
        A matrix is positive definite if it's symmetric and all eigenvalues are positive.
        
        Returns:
            bool: True if the matrix is positive definite, False otherwise.
        """
        if not (self.is_square() and self.is_symmetric()):
            return False
        eigenvalues = self.find_eigenvalues()
        return all(ev > 0 for ev in eigenvalues)

    def norm(self):
        """
        Calculate the Euclidean (Frobenius) norm of the matrix.
        
        Returns:
            float: The Euclidean norm (square root of the sum of squares of all elements).
        """
        return (sum(cell ** 2 for row in self.matrix for cell in row)) ** 0.5
