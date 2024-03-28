from geometry import Geometry, Intersection, epsilon
import glm
import helperclasses as hc
import math


class Sphere(Geometry):
    shadow_epsilon = 10 ** (-4)

    def __init__(self, name: str, gtype: str, materials: list[hc.Material], center: glm.vec3, radius: float, speed):
        super().__init__(name, gtype, materials, speed)
        self.center = center
        self.radius = radius

    def intersect(self, ray: hc.Ray) -> list[Intersection]:
        if self.speed is not None:
            center = self.center + self.speed * self.scene.current_time
        else:
            center = self.center

        a = glm.dot(ray.direction, ray.direction)  # a = d . d
        b = 2 * glm.dot(ray.direction, ray.origin - center)  # b = 2 (d . p - d . c)
        c = glm.dot(ray.origin - center, ray.origin - center) - self.radius ** 2  # c = p . p - 2 p . c + c . c - r^2

        discriminant = b ** 2 - 4 * a * c
        if discriminant < 0:
            return []

        intersections = []
        t1 = (-b - math.sqrt(discriminant)) / (2 * a)
        t2 = (-b + math.sqrt(discriminant)) / (2 * a)
        for t in [t1, t2]:
            if t > 0:
                position = ray.getPoint(t)
                normal = glm.normalize(position - center)
                intersections.append(Intersection(t, normal, position, self.materials[0], self))

        return intersections

    def shadow_intersect(self, ray: hc.Ray, t_max: float) -> bool:
        if self.speed is not None:
            center = self.center + self.speed * self.scene.current_time
        else:
            center = self.center

        a = glm.dot(ray.direction, ray.direction)  # a = d . d
        b = 2 * glm.dot(ray.direction, ray.origin - center)  # b = 2 (d . p - d . c)
        c = glm.dot(ray.origin - center, ray.origin - center) - self.radius ** 2  # c = p . p - 2 p . c + c . c - r^2

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
        if self.speed is not None:
            center = self.center + self.speed * self.scene.current_time
        else:
            center = self.center

        return glm.length(point - center) < self.radius

    def __repr__(self):
        return f"Sphere({self.name}, center: {self.center}, radius: {self.radius})"


class Plane(Geometry):
    def __init__(self, name: str, gtype: str, materials: list[hc.Material], point: glm.vec3, normal: glm.vec3, speed):
        super().__init__(name, gtype, materials, speed)
        self.point = point
        self.normal = normal

    def intersect(self, ray: hc.Ray) -> list[Intersection]:
        if self.speed is not None:
            point = self.point + self.speed * self.scene.current_time
        else:
            point = self.point

        denom = glm.dot(ray.direction, self.normal)
        if abs(denom) > epsilon:
            t = glm.dot(point - ray.origin, self.normal) / denom
            if t >= 0:
                position = ray.getPoint(t)

                if len(self.materials) == 1:
                    mat = self.materials[0]
                else:
                    dx = math.floor(position.x - point.x)
                    dz = math.floor(position.z - point.z)
                    mat = self.materials[(dx + dz) % 2]

                intersect = Intersection(t, self.normal, position, mat, self)
                return [intersect]
        return []

    def shadow_intersect(self, ray: hc.Ray, t_max: float) -> bool:
        if self.speed is not None:
            point = self.point + self.speed * self.scene.current_time
        else:
            point = self.point

        denom = glm.dot(ray.direction, self.normal)
        if abs(denom) > epsilon:
            t = glm.dot(point - ray.origin, self.normal) / denom
            return self.shadow_epsilon < t < t_max

    def get_material(self, point: glm.vec3) -> hc.Material:
        if self.speed is not None:
            point = self.point + self.speed * self.scene.current_time
        else:
            point = self.point

        if len(self.materials) == 1:
            return self.materials[0]
        dx = math.floor(point.x - point.x)
        dz = math.floor(point.z - point.z)
        return self.materials[(dx + dz) % 2]

    def __repr__(self):
        return f"Plane({self.name}, point: {self.point}, normal: {self.normal})"


class AABB(Geometry):
    def __init__(self, name: str, gtype: str, materials: list[hc.Material], center: glm.vec3, dimension: glm.vec3, speed):
        # dimension holds information for length of each size of the box
        super().__init__(name, gtype, materials, speed)
        halfside = dimension / 2
        self.minpos = center - halfside
        self.maxpos = center + halfside

    def intersect(self, ray: hc.Ray) -> list[Intersection]:
        if self.speed is not None:
            minpos = self.minpos + self.speed * self.scene.current_time
            maxpos = self.maxpos + self.speed * self.scene.current_time
        else:
            minpos = self.minpos
            maxpos = self.maxpos

        if ray.direction.x == 0:
            if not (minpos.x < ray.origin.x < maxpos.x):
                return []
            x_interval = hc.AAInterval(float("-inf"), float("inf"), "x")
        else:
            t_x_1 = (minpos.x - ray.origin.x) / ray.direction.x
            t_x_2 = (maxpos.x - ray.origin.x) / ray.direction.x
            x_interval = hc.AAInterval(t_x_1, t_x_2, "x")

        if ray.direction.y == 0:
            if not (minpos.y < ray.origin.y < maxpos.y):
                return []
            y_interval = hc.AAInterval(float("-inf"), float("inf"), "y")
        else:
            t_y_1 = (minpos.y - ray.origin.y) / ray.direction.y
            t_y_2 = (maxpos.y - ray.origin.y) / ray.direction.y
            y_interval = hc.AAInterval(t_y_1, t_y_2, "y")

        if ray.direction.z == 0:
            if not (minpos.z < ray.origin.z < maxpos.z):
                return []
            z_interval = hc.AAInterval(float("-inf"), float("inf"), "z")
        else:
            t_z_1 = (minpos.z - ray.origin.z) / ray.direction.z
            t_z_2 = (maxpos.z - ray.origin.z) / ray.direction.z
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
        if self.speed is not None:
            minpos = self.minpos + self.speed * self.scene.current_time
            maxpos = self.maxpos + self.speed * self.scene.current_time
        else:
            minpos = self.minpos
            maxpos = self.maxpos

        if ray.direction.x == 0:
            if not (minpos.x < ray.origin.x < maxpos.x):
                return False
            x_interval = hc.AAInterval(float("-inf"), float("inf"), "x")
        else:
            t_x_1 = (minpos.x - ray.origin.x) / ray.direction.x
            t_x_2 = (maxpos.x - ray.origin.x) / ray.direction.x
            x_interval = hc.AAInterval(t_x_1, t_x_2, "x")

        if ray.direction.y == 0:
            if not (minpos.y < ray.origin.y < maxpos.y):
                return False
            y_interval = hc.AAInterval(float("-inf"), float("inf"), "y")
        else:
            t_y_1 = (minpos.y - ray.origin.y) / ray.direction.y
            t_y_2 = (maxpos.y - ray.origin.y) / ray.direction.y
            y_interval = hc.AAInterval(t_y_1, t_y_2, "y")

        if ray.direction.z == 0:
            if not (minpos.z < ray.origin.z < maxpos.z):
                return False
            z_interval = hc.AAInterval(float("-inf"), float("inf"), "z")
        else:
            t_z_1 = (minpos.z - ray.origin.z) / ray.direction.z
            t_z_2 = (maxpos.z - ray.origin.z) / ray.direction.z
            z_interval = hc.AAInterval(t_z_1, t_z_2, "z")

        intersected = (max(x_interval, y_interval, z_interval, key=lambda x: x.start),
                       min(x_interval, y_interval, z_interval, key=lambda x: x.end))

        if intersected[0].start > intersected[1].end:
            return False

        time = intersected[0].start

        return self.shadow_epsilon < time < t_max

    def is_inside(self, point: glm.vec3) -> bool:
        if self.speed is not None:
            minpos = self.minpos + self.speed * self.scene.current_time
            maxpos = self.maxpos + self.speed * self.scene.current_time
        else:
            minpos = self.minpos
            maxpos = self.maxpos

        inside_x = minpos.x < point.x < maxpos.x
        inside_y = minpos.y < point.y < maxpos.y
        inside_z = minpos.z < point.z < maxpos.z
        return inside_x and inside_y and inside_z

    def __repr__(self):
        return f"AABB({self.name}, minpos: {self.minpos}, maxpos: {self.maxpos})"
