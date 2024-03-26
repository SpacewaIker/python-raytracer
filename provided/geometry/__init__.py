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
    shadow_epsilon = 10 ** (-6)

    def __init__(self, name: str, gtype: str, materials: list[hc.Material]):
        self.name = name
        self.gtype = gtype
        self.materials = materials

    def intersect(self, ray: hc.Ray) -> Intersection:
        return None

    def shadow_intersect(self, ray: hc.Ray, t_max: float) -> bool:
        return False

    def __repr__(self):
        return f"Geometry({self.name}, type: {self.gtype})"


class Hierarchy(Geometry):
    def __init__(self, name: str, gtype: str, materials: list[hc.Material], t: glm.vec3, r: glm.vec3, s: glm.vec3):
        super().__init__(name, gtype, materials)
        self.t = t
        self.M = glm.mat4(1.0)
        self.Minv = glm.mat4(1.0)
        self.make_matrices(t, r, s)
        self.children: list[Geometry] = []

    def make_matrices(self, t: glm.vec3, r: glm.vec3, s: glm.vec3):
        self.M = glm.mat4(1.0)
        self.M = glm.translate(self.M, t)
        self.M = glm.rotate(self.M, glm.radians(r.x), glm.vec3(1, 0, 0))
        self.M = glm.rotate(self.M, glm.radians(r.y), glm.vec3(0, 1, 0))
        self.M = glm.rotate(self.M, glm.radians(r.z), glm.vec3(0, 0, 1))
        self.M = glm.scale(self.M, s)
        self.Minv = glm.inverse(self.M)
        self.t = t
        self.r = r
        self.s = s

    def intersect(self, ray: hc.Ray) -> Intersection:
        m_inv_o = glm.vec3(self.Minv * glm.vec4(ray.origin, 1))
        m_inv_d = glm.vec3(self.Minv * glm.vec4(ray.direction, 0))
        m_ray = hc.Ray(m_inv_o, m_inv_d)

        intersections = []
        for child in self.children:
            intersect = child.intersect(m_ray)
            if intersect is not None and intersect.time < float('inf'):
                intersections.append(intersect)

        if len(intersections) == 0:
            return None

        first_intersect = min(intersections, key=lambda x: x.time)

        if first_intersect.mat is None:
            first_intersect.mat = self.materials[0]

        first_intersect.position = (self.M * glm.vec4(first_intersect.position, 1)).xyz
        first_intersect.normal = glm.normalize(glm.transpose(self.Minv) * glm.vec4(first_intersect.normal, 0)).xyz

        return first_intersect

    def shadow_intersect(self, ray: hc.Ray, t_max: float) -> bool:
        m_inv_o = glm.vec3(self.Minv * glm.vec4(ray.origin, 1))
        m_inv_d = glm.vec3(self.Minv * glm.vec4(ray.direction, 0))
        m_ray = hc.Ray(m_inv_o, m_inv_d)

        for child in self.children:
            if child.shadow_intersect(m_ray, t_max):
                return True

        return False

    def __repr__(self):
        return f"Hierarchy({self.name}, t: {self.t}, r: {self.r}, s: {self.s}, children: {len(self.children)})"
