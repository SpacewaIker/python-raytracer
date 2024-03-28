from geometry import Geometry, Intersection, epsilon
import glm
import helperclasses as hc
import math


class Sphere(Geometry):
    shadow_epsilon = 10 ** (-4)

    def __init__(self, name: str, gtype: str, materials: list[hc.Material], center: glm.vec3, radius: float):
        super().__init__(name, gtype, materials)
        self.center = center
        self.radius = radius

    def intersect(self, ray: hc.Ray) -> list[Intersection]:
        a = glm.dot(ray.direction, ray.direction)  # a = d . d
        b = 2 * glm.dot(ray.direction, ray.origin -
                        self.center)  # b = 2 (d . p - d . c)
        c = glm.dot(ray.origin - self.center, ray.origin - self.center) - \
            self.radius ** 2  # c = p . p - 2 p . c + c . c - r^2

        discriminant = b ** 2 - 4 * a * c
        if discriminant < 0:
            return []

        intersections = []
        t1 = (-b - math.sqrt(discriminant)) / (2 * a)
        t2 = (-b + math.sqrt(discriminant)) / (2 * a)
        for t in [t1, t2]:
            if t > 0:
                position = ray.getPoint(t)
                normal = glm.normalize(position - self.center)
                intersections.append(Intersection(t, normal, position, self.materials[0], self))

        return intersections

    def shadow_intersect(self, ray: hc.Ray, t_max: float) -> bool:
        a = glm.dot(ray.direction, ray.direction)  # a = d . d
        b = 2 * glm.dot(ray.direction, ray.origin -
                        self.center)  # b = 2 (d . p - d . c)
        c = glm.dot(ray.origin - self.center, ray.origin - self.center) - \
            self.radius ** 2  # c = p . p - 2 p . c + c . c - r^2

        discriminant = b ** 2 - 4 * a * c
        if discriminant < 0:
            return False

        t = (-b - math.sqrt(discriminant)) / (2 * a)

        if self.shadow_epsilon < t < t_max:
            return True

        t = (-b + math.sqrt(discriminant)) / (2 * a)

        if self.shadow_epsilon < t < t_max:
            return True

        return False

    def is_inside(self, point: glm.vec3) -> bool:
        return glm.length(point - self.center) < self.radius

    def __repr__(self):
        return f"Sphere({self.name}, center: {self.center}, radius: {self.radius})"


class Plane(Geometry):
    def __init__(self, name: str, gtype: str, materials: list[hc.Material], point: glm.vec3, normal: glm.vec3):
        super().__init__(name, gtype, materials)
        self.point = point
        self.normal = normal

    def intersect(self, ray: hc.Ray) -> list[Intersection]:
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

                intersect = Intersection(t, self.normal, position, mat, self)
                return [intersect]
        return []

    def shadow_intersect(self, ray: hc.Ray, t_max: float) -> bool:
        denom = glm.dot(ray.direction, self.normal)
        if abs(denom) > epsilon:
            t = glm.dot(self.point - ray.origin, self.normal) / denom
            return self.shadow_epsilon < t < t_max

    def get_material(self, point: glm.vec3) -> hc.Material:
        if len(self.materials) == 1:
            return self.materials[0]
        dx = math.floor(point.x - self.point.x)
        dz = math.floor(point.z - self.point.z)
        return self.materials[(dx + dz) % 2]

    def __repr__(self):
        return f"Plane({self.name}, point: {self.point}, normal: {self.normal})"


class AABB(Geometry):
    def __init__(self, name: str, gtype: str, materials: list[hc.Material], center: glm.vec3, dimension: glm.vec3):
        # dimension holds information for length of each size of the box
        super().__init__(name, gtype, materials)
        halfside = dimension / 2
        self.minpos = center - halfside
        self.maxpos = center + halfside

    def intersect(self, ray: hc.Ray) -> list[Intersection]:
        if ray.direction.x == 0:
            if not (self.minpos.x < ray.origin.x < self.maxpos.x):
                return []
            x_interval = hc.AAInterval(float("-inf"), float("inf"), "x")
        else:
            t_x_1 = (self.minpos.x - ray.origin.x) / ray.direction.x
            t_x_2 = (self.maxpos.x - ray.origin.x) / ray.direction.x
            x_interval = hc.AAInterval(t_x_1, t_x_2, "x")

        if ray.direction.y == 0:
            if not (self.minpos.y < ray.origin.y < self.maxpos.y):
                return []
            y_interval = hc.AAInterval(float("-inf"), float("inf"), "y")
        else:
            t_y_1 = (self.minpos.y - ray.origin.y) / ray.direction.y
            t_y_2 = (self.maxpos.y - ray.origin.y) / ray.direction.y
            y_interval = hc.AAInterval(t_y_1, t_y_2, "y")

        if ray.direction.z == 0:
            if not (self.minpos.z < ray.origin.z < self.maxpos.z):
                return []
            z_interval = hc.AAInterval(float("-inf"), float("inf"), "z")
        else:
            t_z_1 = (self.minpos.z - ray.origin.z) / ray.direction.z
            t_z_2 = (self.maxpos.z - ray.origin.z) / ray.direction.z
            z_interval = hc.AAInterval(t_z_1, t_z_2, "z")

        intersected = (max(x_interval, y_interval, z_interval, key=lambda x: x.start),
                       min(x_interval, y_interval, z_interval, key=lambda x: x.end))

        if intersected[0].start > intersected[1].end or intersected[0].start < 0:
            return []

        intersections = []
        for t in [intersected[0].start, intersected[1].end]:
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

            position = ray.getPoint(t)
            mat = self.materials[0]

            intersections.append(Intersection(t, normal, position, mat, self))

        return intersections

    def shadow_intersect(self, ray: hc.Ray, t_max: float) -> bool:
        if ray.direction.x == 0:
            if not (self.minpos.x < ray.origin.x < self.maxpos.x):
                return False
            x_interval = hc.AAInterval(float("-inf"), float("inf"), "x")
        else:
            t_x_1 = (self.minpos.x - ray.origin.x) / ray.direction.x
            t_x_2 = (self.maxpos.x - ray.origin.x) / ray.direction.x
            x_interval = hc.AAInterval(t_x_1, t_x_2, "x")

        if ray.direction.y == 0:
            if not (self.minpos.y < ray.origin.y < self.maxpos.y):
                return False
            y_interval = hc.AAInterval(float("-inf"), float("inf"), "y")
        else:
            t_y_1 = (self.minpos.y - ray.origin.y) / ray.direction.y
            t_y_2 = (self.maxpos.y - ray.origin.y) / ray.direction.y
            y_interval = hc.AAInterval(t_y_1, t_y_2, "y")

        if ray.direction.z == 0:
            if not (self.minpos.z < ray.origin.z < self.maxpos.z):
                return False
            z_interval = hc.AAInterval(float("-inf"), float("inf"), "z")
        else:
            t_z_1 = (self.minpos.z - ray.origin.z) / ray.direction.z
            t_z_2 = (self.maxpos.z - ray.origin.z) / ray.direction.z
            z_interval = hc.AAInterval(t_z_1, t_z_2, "z")

        intersected = (max(x_interval, y_interval, z_interval, key=lambda x: x.start),
                       min(x_interval, y_interval, z_interval, key=lambda x: x.end))

        if intersected[0].start > intersected[1].end:
            return False

        time = intersected[0].start

        return self.shadow_epsilon < time < t_max

    def is_inside(self, point: glm.vec3) -> bool:
        inside_x = self.minpos.x < point.x < self.maxpos.x
        inside_y = self.minpos.y < point.y < self.maxpos.y
        inside_z = self.minpos.z < point.z < self.maxpos.z
        return inside_x and inside_y and inside_z

    def __repr__(self):
        return f"AABB({self.name}, minpos: {self.minpos}, maxpos: {self.maxpos})"
