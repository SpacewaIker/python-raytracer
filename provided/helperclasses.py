import taichi.math as tm
import taichi as ti

# Ported from C++ by Melissa Katz
# Adapted from code by Loïc Nassif and Paul Kry

vec3 = ti.math.vec3


@ti.dataclass
class Ray:
    origin: vec3
    direction: vec3

    @ti.func
    def __init__(self, o: vec3, d: vec3):
        self.origin = o
        self.direction = d

    @ti.func
    def getDistance(self, point: vec3):
        return tm.length(point - self.origin)

    @ti.func
    def getPoint(self, t: float):
        return self.origin + self.direction * t


@ti.dataclass
class Material:
    specular: vec3
    diffuse: vec3
    hardness: float
    ID: int

    def __init__(self, specular: vec3, diffuse: vec3, hardness: float, ID: int):
        self.specular = specular
        self.diffuse = diffuse
        self.hardness = hardness
        self.ID = ID

    @staticmethod
    def default():
        specular = diffuse = vec3(0, 0, 0)
        hardness = ID = -1
        return Material(specular, diffuse, hardness, ID)


@ti.dataclass
class Light:
    ltype: int
    colour: vec3
    vector: vec3
    power: float

    def __init__(self, ltype: int, colour: vec3, vector: vec3, power: float):
        self.ltype = ltype
        self.colour = colour
        self.vector = vector
        self.power = power


@ti.dataclass
class Intersection:
    time: float
    normal: vec3
    position: vec3
    mat: Material

    def __init__(self, time: float, normal: vec3, position: vec3, material: Material):
        self.time = time
        self.normal = normal
        self.position = position
        self.mat = material

    @staticmethod
    def default():
        time = float("inf")
        normal = vec3(0, 0, 0)
        position = vec3(0, 0, 0)
        mat = Material.default()
        return Intersection(time, normal, position, mat)
