import math
import helperclasses as hc
import glm
import igl

# Ported from C++ by Melissa Katz
# Adapted from code by Loïc Nassif and Paul Kry

epsilon = 10 ** (-4)


class Geometry:
    def __init__(self, name: str, gtype: str, materials: list[hc.Material]):
        self.name = name
        self.gtype = gtype
        self.materials = materials

    def intersect(self, ray: hc.Ray, intersect: hc.Intersection):
        return intersect


class Sphere(Geometry):
    def __init__(self, name: str, gtype: str, materials: list[hc.Material], center: glm.vec3, radius: float):
        super().__init__(name, gtype, materials)
        self.center = center
        self.radius = radius

    def intersect(self, ray: hc.Ray, intersect: hc.Intersection):
        a = glm.dot(ray.direction, ray.direction)  # a = d . d
        b = 2 * glm.dot(ray.direction, ray.origin -
                        self.center)  # b = 2 (d . p - d . c)
        c = glm.dot(ray.origin - self.center, ray.origin - self.center) - \
            self.radius ** 2  # c = p . p - 2 p . c + c . c - r^2

        discriminant = b ** 2 - 4 * a * c
        if discriminant < 0:
            return None
        t = (-b - math.sqrt(discriminant)) / (2 * a)
        position = ray.getPoint(t)
        normal = glm.normalize(position - self.center)

        intersect = hc.Intersection(t, normal, position, self.materials[0])

        return intersect


class Plane(Geometry):
    def __init__(self, name: str, gtype: str, materials: list[hc.Material], point: glm.vec3, normal: glm.vec3):
        super().__init__(name, gtype, materials)
        self.point = point
        self.normal = normal

    def intersect(self, ray: hc.Ray, intersect: hc.Intersection):
        denom = glm.dot(ray.direction, self.normal)
        if abs(denom) > epsilon:
            t = glm.dot(self.point - ray.origin, self.normal) / denom
            if t >= 0:
                position = ray.getPoint(t)

                if len(self.materials) == 1:
                    mat = self.materials[0]
                else:
                    dx = math.floor(position.x - self.point.x)
                    dz = math.floor(position.z - self.point.z)
                    mat = self.materials[(dx + dz) % 2]

                intersect = hc.Intersection(t, self.normal, position, mat)
                return intersect
        return None


class AABB(Geometry):
    def __init__(self, name: str, gtype: str, materials: list[hc.Material], center: glm.vec3, dimension: glm.vec3):
        # dimension holds information for length of each size of the box
        super().__init__(name, gtype, materials)
        halfside = dimension / 2
        self.minpos = center - halfside
        self.maxpos = center + halfside

    def intersect(self, ray: hc.Ray, intersect: hc.Intersection):
        if ray.direction.x == 0:
            if not (self.minpos.x < ray.origin.x < self.maxpos.x):
                return None
            x_interval = hc.AAInterval(float("-inf"), float("inf"), "x")
        else:
            t_x_1 = (self.minpos.x - ray.origin.x) / ray.direction.x
            t_x_2 = (self.maxpos.x - ray.origin.x) / ray.direction.x
            x_interval = hc.AAInterval(t_x_1, t_x_2, "x")

        if ray.direction.y == 0:
            if not (self.minpos.y < ray.origin.y < self.maxpos.y):
                return None
            y_interval = hc.AAInterval(float("-inf"), float("inf"), "y")
        else:
            t_y_1 = (self.minpos.y - ray.origin.y) / ray.direction.y
            t_y_2 = (self.maxpos.y - ray.origin.y) / ray.direction.y
            y_interval = hc.AAInterval(t_y_1, t_y_2, "y")

        if ray.direction.z == 0:
            if not (self.minpos.z < ray.origin.z < self.maxpos.z):
                return None
            z_interval = hc.AAInterval(float("-inf"), float("inf"), "z")
        else:
            t_z_1 = (self.minpos.z - ray.origin.z) / ray.direction.z
            t_z_2 = (self.maxpos.z - ray.origin.z) / ray.direction.z
            z_interval = hc.AAInterval(t_z_1, t_z_2, "z")

        intersected = (max(x_interval, y_interval, z_interval, key=lambda x: x.start),
                       min(x_interval, y_interval, z_interval, key=lambda x: x.end))

        if intersected[0].start > intersected[1].end:
            return None

        time = intersected[0].start

        if intersected[0].label == "x" and ray.direction.x < 0:
            normal = glm.vec3(1, 0, 0)
        elif intersected[0].label == "x" and ray.direction.x > 0:
            normal = glm.vec3(-1, 0, 0)
        elif intersected[0].label == "y" and ray.direction.y < 0:
            normal = glm.vec3(0, 1, 0)
        elif intersected[0].label == "y" and ray.direction.y > 0:
            normal = glm.vec3(0, -1, 0)
        elif intersected[0].label == "z" and ray.direction.z < 0:
            normal = glm.vec3(0, 0, 1)
        elif intersected[0].label == "z" and ray.direction.z > 0:
            normal = glm.vec3(0, 0, -1)

        position = ray.getPoint(time)
        mat = self.materials[0]

        intersect = hc.Intersection(time, normal, position, mat)

        return intersect


class Mesh(Geometry):
    def __init__(self, name: str, gtype: str, materials: list[hc.Material], translate: glm.vec3, scale: float,
                 filepath: str):
        super().__init__(name, gtype, materials)
        verts, _, norms, self.faces, _, _ = igl.read_obj(filepath)
        self.verts = []
        self.norms = []
        for v in verts:
            self.verts.append((glm.vec3(v[0], v[1], v[2]) + translate) * scale)
        for n in norms:
            self.norms.append(glm.vec3(n[0], n[1], n[2]))

    def intersect(self, ray: hc.Ray, intersect: hc.Intersection):
        pass
        # TODO: Create intersect code for Mesh


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

    def intersect(self, ray: hc.Ray, intersect: hc.Intersection):
        pass
        # TODO: Create intersect code for Hierarchy
