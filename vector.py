class Vector:
    """
    A class representing a mathematical vector with various operations.
    
    Attributes:
        vector (list): A list containing the vector components.
        size (int): The dimension of the vector.
    """
    def __init__(self, vector):
        """
        Initialize a vector object from a list of components.
        
        Parameters:
            vector (list): A list of the vector components.
        """
        self.vector = vector
        self.size = len(vector)

    def show_vector(self):
        """
        Display the vector components.
        Prints each element of the vector followed by a space.
        """
        for i in range(self.size):
            print(self.vector[i], end=' ')
        print()
    
    def add(self, other):
        """
        Add another vector to this vector.
        
        Parameters:
            other (Vector): The vector to add to this vector.
            
        Returns: 
            Vector: The resulting vector after addition.

        Raises:
            ValueError: If the vectors have different dimensions
        """
        if self.size != other.size:
            raise ValueError("Vectors must have same dimensions")
        ans = [0 for _ in range(self.size)]
        for i in range(self.size):
            ans[i] = self.vector[i] + other.vector[i]
        return Vector(ans)

    def subtract(self, other):
        """
        Subtract another vector from this vector.
        
        Parameters:
            other (Vector): The vector to subtract.
            
        Returns:
            Vector: The resulting vector after subtraction.
            
        Raises:
            ValueError: If the vectors have different dimensions
        """
        if self.size != other.size:
            raise ValueError("Vectors must have same dimensions")
        ans = [0 for _ in range(self.size)]
        for i in range(self.size):
            ans[i] = self.vector[i] - other.vector[i]
        return Vector(ans)

    def dot_product(self, other) -> float:
        """
        Calculate the dot product of this vector with another vector.
        
        Parameters:
            other (Vector): The vector to calculate the dot product with.
            
        Returns:
            float: The dot product of the two vectors.
            
        Raises:
            ValueError: If the vectors have different dimensions.
        """
        if self.size != other.size:
            raise ValueError('Vectors must have same dimensions')
        dotProduct = 0
        for i in range(self.size):
            dotProduct += self.vector[i] * other.vector[i]
        return float(dotProduct)

    def scalar_multiply(self, k: float|int):
        """
        Multiply the vector by a scalar.
        
        Parameters:
            k (float | int): The scalar value to multiply by.
            
        Returns:
            Vector: The resulting vector after scalar multiplication.
        """
        return Vector([k * x for x in self.vector])

    def cross_product(self, other):
        """
        Calculate the cross product of this vector with another vector.
        
        For 2D vectors: Returns a scalar (the z-component of the cross product).
        For 3D vectors: Returns a Vector representing the full cross product.
        
        Parameters:
            other (Vector): The vector to calculate the cross product with.
            
        Returns:
            Vector|float: Cross product result (Vector for 3D, scalar for 2D).
            
        Raises:
            ValueError: If the vectors have different dimensions or if dimensions are not 2 or 3.
        """
        if self.size != other.size:
            raise ValueError("Vectors must have the same dimensions")
        if self.size == 2:
            return self.vector[0] * other.vector[1] - self.vector[1] * other.vector[0]
        elif self.size == 3:
            result = [
                self.vector[1] * other.vector[2] - self.vector[2] * other.vector[1],
                self.vector[2] * other.vector[0] - self.vector[0] * other.vector[2],
                self.vector[0] * other.vector[1] - self.vector[1] * other.vector[0]
            ]
            return Vector(result)
        else:
            raise ValueError("Cross product is defined only for 2D and 3D vectors")

    def norm(self, p=2) -> float:
        """
        Calculate the p-norm of the vector.
        
        Parameters:
            p (int, default=2): The order of the norm. p=2 gives the Euclidean (L2) norm.
            
        Returns:
            float: The p-norm of the vector.
        """
        return (sum(abs(x)**p for x in self.vector))**(1/p)
