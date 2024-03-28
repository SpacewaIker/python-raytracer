import helperclasses as hc
import glm

# Ported from C++ by Melissa Katz
# Adapted from code by Loïc Nassif and Paul Kry

epsilon = 10 ** (-4)


class Intersection:

    def __init__(self, time: float, normal: glm.vec3, position: glm.vec3, material: hc.Material, geometry: "Geometry" = None):
        self.time = time
        self.normal = normal
        self.position = position
        self.mat = material
        self.geometry = geometry

    @staticmethod
    def default():
        time = float("inf")
        normal = glm.vec3(0, 0, 0)
        position = glm.vec3(0, 0, 0)
        mat = hc.Material.default()
        return Intersection(time, normal, position, mat)

    def __repr__(self):
        if self.geometry is None:
            return f"Intersection(time: {self.time}, position: {self.position}, normal: {self.normal})"
        return f"Intersection(time: {self.time}, position: {self.position}, normal: {self.normal}, geometry: {self.geometry.name})"


class Geometry:
    shadow_epsilon = 10 ** (-5)

    def __init__(self, name: str, gtype: str, materials: list[hc.Material], speed):
        self.name = name
        self.gtype = gtype
        self.materials = materials
        self.speed = speed

    def intersect(self, ray: hc.Ray) -> list[Intersection]:
        return []

    def shadow_intersect(self, ray: hc.Ray, t_max: float) -> bool:
        return False

    def is_inside(self, point: glm.vec3) -> bool:
        return False

    def get_material(self, point: glm.vec3) -> hc.Material:
        return self.materials[0]

    def set_scene(self, scene):
        self.scene = scene

    def __repr__(self):
        return f"Geometry({self.name}, type: {self.gtype})"
