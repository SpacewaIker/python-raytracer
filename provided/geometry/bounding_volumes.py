# Thibaut Baguette
# 261001513
# COMP557/ECSE532 Assignment 4
# 2024-04-08

import glm
from geometry import Geometry
import math
import helperclasses as hc


class BoundingSphere():
    def __init__(self, center: glm.vec3, radius: float, geometry: Geometry):
        self.center = center
        self.radius = radius
        self.geometry = geometry

    def intersect(self, ray: hc.Ray) -> bool:
        a = glm.dot(ray.direction, ray.direction)  # a = d . d
        b = 2 * glm.dot(ray.direction, ray.origin -
                        self.center)  # b = 2 (d . p - d . c)
        c = glm.dot(ray.origin - self.center, ray.origin - self.center) - \
            self.radius ** 2  # c = p . p - 2 p . c + c . c - r^2

        discriminant = b ** 2 - 4 * a * c
        if discriminant < 0:
            return False

        t = (-b - math.sqrt(discriminant)) / (2 * a)
        if t > 0:
            return True

        t = (-b + math.sqrt(discriminant)) / (2 * a)
        if t > 0:
            return True

        return False

    def __repr__(self):
        return f"BoundingSphere({self.geometry.name})"


class BoundingAABB():
    def __init__(self, minpos: glm.vec3, maxpos: glm.vec3, geometry: Geometry):
        # dimension holds information for length of each size of the box
        self.minpos = minpos
        self.maxpos = maxpos

    def intersect(self, ray: hc.Ray) -> bool:
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

        if intersected[0].start > intersected[1].end or intersected[0].start < 0:
            return False

        return True

    def __repr__(self):
        return f"BoundingAABB({self.name}, minpos: {self.minpos}, maxpos: {self.maxpos})"
