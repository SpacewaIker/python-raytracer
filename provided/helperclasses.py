import glm

# Ported from C++ by Melissa Katz
# Adapted from code by Loïc Nassif and Paul Kry


class Ray:
    def __init__(self, o: glm.vec3, d: glm.vec3):
        self.origin = o
        self.direction = d

    def getDistance(self, point: glm.vec3):
        return glm.length(point - self.origin)

    def getPoint(self, t: float):
        return self.origin + self.direction * t


class Material:
    def __init__(self, name: str, specular: glm.vec3, diffuse: glm.vec3, hardness: float, ID: int):
        self.name = name
        self.specular = specular
        self.diffuse = diffuse
        self.hardness = hardness
        self.ID = ID

    @staticmethod
    def default():
        name = "default"
        specular = diffuse = glm.vec3(0, 0, 0)
        hardness = ID = -1
        return Material(name, specular, diffuse, hardness, ID)


class Light:
    def __init__(self, ltype: str, name: str, colour: glm.vec3, vector: glm.vec3, power: float):
        self.type = ltype
        self.name = name
        self.colour = colour
        self.vector = vector
        self.power = power


class Intersection:

    def __init__(self, time: float, normal: glm.vec3, position: glm.vec3, material: Material):
        self.time = time
        self.normal = normal
        self.position = position
        self.mat = material

    @staticmethod
    def default():
        time = float("inf")
        normal = glm.vec3(0, 0, 0)
        position = glm.vec3(0, 0, 0)
        mat = Material.default()
        return Intersection(time, normal, position, mat)


class AAInterval:
    def __init__(self, t1: float, t2: float, label: str = None):
        self.start = min(t1, t2)
        self.end = max(t1, t2)
        self.label = label
