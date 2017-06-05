import math

class gplot_vector3:
    def __init__(self, x, y, z, dx, dy, dz):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.dx = float(dx)
        self.dy = float(dy)
        self.dz = float(dz)

    def __str__(self):
        return str(self.x) + " " + str(self.y) + " " + str(self.z) + " " + str(self.dx) + " " + str(self.dy) + " " +  str(self.dz)

    def __abs__(self):
        return math.sqrt(self.dx**2 + self.dy**2 + self.dz**2)
