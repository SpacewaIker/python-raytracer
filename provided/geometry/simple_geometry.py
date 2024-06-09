# Thibaut Baguette
# 261001513
# COMP557/ECSE532 Assignment 4
# 2024-04-08

from geometry import Geometry, Intersection, epsilon
import glm
import helperclasses as hc
import math


class Sphere(Geometry):
    shadow_epsilon = 10 ** (-3)

    def __init__(self, name: str, gtype: str, materials: list[hc.Material], center: glm.vec3, radius: float, speed):
        super().__init__(name, gtype, materials, speed)
        self.center = center
        self.radius = radius

    def intersect(self, ray: hc.Ray) -> list[Intersection]:
        if self.speed is not None:
            center = self.center + self.speed * self.scene.current_time
        else:
            center = self.center

        if abs(ray.direction.x) < 0.01 and abs(ray.direction.y) < 0.01 and abs(ray.direction.z) < 0.01:
            print("Ray direction is zero ", ray.direction)

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
        self.texture = None

        if self.normal == glm.vec3(0, 1, 0) or self.normal == glm.vec3(0, -1, 0) or self.normal == glm.vec3(0, 0, 1):
            self.width_axis = glm.vec3(1, 0, 0)
        elif self.normal == glm.vec3(0, 0, -1):
            self.width_axis = glm.vec3(-1, 0, 0)
        elif self.normal == glm.vec3(1, 0, 0):
            self.width_axis = glm.vec3(0, 0, -1)
        elif self.normal == glm.vec3(-1, 0, 0):
            self.width_axis = glm.vec3(0, 0, 1)
        else:
            self.width_axis = glm.normalize(glm.cross(self.normal, glm.vec3(0, 0, 1)))
        self.height_axis = glm.normalize(glm.cross(self.width_axis, self.normal))

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
                mat = self.get_material(position)

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
            position = self.point + self.speed * self.scene.current_time
        else:
            position = self.point

        if len(self.materials) == 1:
            return self.materials[0]

        point = point - glm.dot(point - position, self.normal) * self.normal
        x = glm.dot(point - position, self.width_axis)
        z = glm.dot(point - position, self.height_axis)

        dx = math.floor(position.x - x)
        dz = math.floor(position.z - z)
        return self.materials[(dx + dz) % 2]

    def get_diffuse(self, point: glm.vec3) -> glm.vec3:
        if self.texture is not None:
            if self.speed is not None:
                position = self.point + self.speed * self.scene.current_time
            else:
                position = self.point

            try:
                scale = self.texture_scale
            except AttributeError:
                scale = 1.0

            # project point onto plane using plane normal
            point = point - glm.dot(point - self.point, self.normal) * self.normal

            i = int((glm.dot(point - position, self.width_axis) * 1000.0 / scale) % self.texture.width)
            j = int((glm.dot(point - position, self.height_axis) * 1000.0 / scale) % self.texture.height)

            pixel = self.texture.getpixel((i, j))
            colour = glm.vec3(pixel[0] / 255, pixel[1] / 255, pixel[2] / 255)

            return colour
        else:
            return self.get_material(point).diffuse

    def __repr__(self):
        return f"Plane({self.name}, point: {self.point}, normal: {self.normal})"


class AABB(Geometry):
    def __init__(self, name: str, gtype: str, materials: list[hc.Material], center: glm.vec3, dimension: glm.vec3, speed):
        # dimension holds information for length of each size of the box
        super().__init__(name, gtype, materials, speed)
        halfside = dimension / 2
        self.minpos = center - halfside
        self.maxpos = center + halfside
        self.texture = None

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

    def get_diffuse(self, point: glm.vec3) -> glm.vec3:
        if self.texture is not None:
            if self.speed is not None:
                minpos = self.minpos + self.speed * self.scene.current_time
                maxpos = self.maxpos + self.speed * self.scene.current_time
            else:
                minpos = self.minpos
                maxpos = self.maxpos

            x = (point.x - minpos.x) / (maxpos.x - minpos.x)
            y = (point.y - minpos.y) / (maxpos.y - minpos.y)
            z = (point.z - minpos.z) / (maxpos.z - minpos.z)

            if abs(point.x - minpos.x) < epsilon:  # negative x
                i = z * self.texture.width
                j = (1 - y) * self.texture.height
            elif abs(point.x - maxpos.x) < epsilon:  # positive x
                i = (1 - z) * self.texture.width
                j = (1 - y) * self.texture.height
            elif abs(point.y - minpos.y) < epsilon:  # negative y
                i = x * self.texture.width
                j = (1 - z) * self.texture.height
            elif abs(point.y - maxpos.y) < epsilon:  # positive y
                i = x * self.texture.width
                j = z * self.texture.height
            elif abs(point.z - minpos.z) < epsilon:  # negative z
                i = (1 - x) * self.texture.width
                j = (1 - y) * self.texture.height
            elif abs(point.z - maxpos.z) < epsilon:  # positive z
                i = x * self.texture.width
                j = (1 - y) * self.texture.height
            else:
                i = 0
                j = 0

            i = min(max(0, i), self.texture.width - 1)
            j = min(max(0, j), self.texture.height - 1)

            pixel = self.texture.getpixel((i, j))
            colour = glm.vec3(pixel[0] / 255, pixel[1] / 255, pixel[2] / 255)

            return colour
        else:
            return self.materials[0].diffuse
