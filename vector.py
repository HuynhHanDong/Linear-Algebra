class Vector:
    def __init__(self, vector):
        self.vector = vector
        self.size = len(vector)

    def show_vector(self):
        for i in range(self.size):
            print(self.vector[i], end=' ')
        print()
    
    def add(self, other):
        if self.size != other.size:
            raise ValueError("Vectors must have same dimensions")
        ans = [0 for _ in range(self.size)]
        for i in range(self.size):
            ans[i] = self.vector[i] + other.vector[i]
        return Vector(ans)

    def subtract(self, other):
        if self.size != other.size:
            raise ValueError("Vectors must have same dimensions")
        ans = [0 for _ in range(self.size)]
        for i in range(self.size):
            ans[i] = self.vector[i] - other.vector[i]
        return Vector(ans)

    def dot_product(self, other):
        if self.size != other.size:
            print('Vectors must have same dimensions')
            return
        dotProduct = 0
        for i in range(self.size):
            dotProduct += self.vector[i] * other.vector[i]
        return dotProduct

    def scalar_multiply(self, k):
        return Vector([k * x for x in self.vector])

    def cross_product(self, other):
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

    def norm(self, p=2):
        return (sum(abs(x)**p for x in self.vector))**(1/p)
