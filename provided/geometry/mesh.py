from geometry import Geometry, Intersection, epsilon
import helperclasses as hc
import glm
import igl
import math
from geometry.bounding_volumes import BoundingAABB, BoundingSphere
from collections import defaultdict
import numpy as np


class Mesh(Geometry):
    def __init__(self, name: str, gtype: str, materials: list[hc.Material], translate: glm.vec3, scale: float,
                 filepath: str, flat_shaded: bool = False):
        super().__init__(name, gtype, materials)
        verts, _, norms, self.faces, _, _ = igl.read_obj(filepath)
        self.verts = []
        self.norms = []
        for v in verts:
            self.verts.append((glm.vec3(v[0], v[1], v[2]) + translate) * scale)
        for n in norms:
            self.norms.append(glm.vec3(n[0], n[1], n[2]))
        if not flat_shaded:
            self._compute_normals()
        self.flat_shaded = flat_shaded

        max_x = max([v.x for v in self.verts])
        min_x = min([v.x for v in self.verts])
        max_y = max([v.y for v in self.verts])
        min_y = min([v.y for v in self.verts])
        max_z = max([v.z for v in self.verts])
        min_z = min([v.z for v in self.verts])

        avg_x = (max_x + min_x) / 2
        avg_y = (max_y + min_y) / 2
        avg_z = (max_z + min_z) / 2
        center = glm.vec3(avg_x, avg_y, avg_z)

        max_dist = max([glm.length(v - center) for v in self.verts])

        aabb_volume = (max_x - min_x) * (max_y - min_y) * (max_z - min_z)
        sphere_volume = 4 / 3 * math.pi * max_dist ** 3

        if aabb_volume < sphere_volume:
            self.bounding_volume = BoundingAABB(glm.vec3(min_x, min_y, min_z), glm.vec3(max_x, max_y, max_z), None)
        else:
            self.bounding_volume = BoundingSphere(center, max_dist, None)

    def _compute_normals(self):
        normals = defaultdict(glm.vec3)
        for face in self.faces:
            v0 = self.verts[face[0]]
            v1 = self.verts[face[1]]
            v2 = self.verts[face[2]]

            e1 = v1 - v0
            e2 = v2 - v0
            normal = glm.normalize(glm.cross(e1, e2))
            area = glm.length(glm.cross(e1, e2)) / 2
            weighted_normal = normal * area

            normals[face[0]] += weighted_normal
            normals[face[1]] += weighted_normal
            normals[face[2]] += weighted_normal

        self.norms = [glm.normalize(normals[i]) for i in range(len(self.verts))]

    def intersect(self, ray: hc.Ray) -> Intersection:
        if not self.bounding_volume.intersect(ray):
            return None

        intersections = []
        for face in self.faces:
            v0 = self.verts[face[0]]
            v1 = self.verts[face[1]]
            v2 = self.verts[face[2]]

            e1 = v1 - v0
            e2 = v2 - v0
            normal = glm.normalize(glm.cross(e1, e2))

            # check plane intersection
            denom = glm.dot(ray.direction, normal)
            if abs(denom) < epsilon:
                continue

            time = glm.dot(v0 - ray.origin, normal) / denom

            if time < 0:
                continue

            # check if intersection point is inside triangle
            point = ray.getPoint(time)
            b0 = glm.dot(glm.cross(v1 - v0, point - v0), normal)
            b1 = glm.dot(glm.cross(v2 - v1, point - v1), normal)
            b2 = glm.dot(glm.cross(v0 - v2, point - v2), normal)

            if b0 >= 0 and b1 >= 0 and b2 >= 0:
                if not self.flat_shaded:
                    p = np.array(point).reshape(1, 3)
                    v0 = np.array(v0).reshape(1, 3)
                    v1 = np.array(v1).reshape(1, 3)
                    v2 = np.array(v2).reshape(1, 3)
                    bar = igl.barycentric_coordinates_tri(p, v0, v1, v2)
                    b0 = bar[0]
                    b1 = bar[1]
                    b2 = bar[2]

                    normal = glm.normalize(b0 * self.norms[face[0]] + b1 * self.norms[face[1]] + b2 * self.norms[face[2]])

                mat = self.materials[0]
                intersect = Intersection(time, normal, point, mat, self)
                intersections.append(intersect)

        if len(intersections) == 0:
            return None

        return min(intersections, key=lambda x: x.time)

    def shadow_intersect(self, ray: hc.Ray, t_max: float) -> bool:
        if not self.bounding_volume.intersect(ray):
            return False

        for face in self.faces:
            v0 = self.verts[face[0]]
            v1 = self.verts[face[1]]
            v2 = self.verts[face[2]]

            e1 = v1 - v0
            e2 = v2 - v0
            normal = glm.cross(e1, e2)

            # check plane intersection
            denom = glm.dot(ray.direction, normal)
            if abs(denom) < epsilon:
                continue

            time = glm.dot(v0 - ray.origin, normal) / denom

            if time < self.shadow_epsilon:
                continue

            # check if intersection point is inside triangle
            point = ray.getPoint(time)
            c0 = glm.cross(v1 - v0, point - v0)
            c1 = glm.cross(v2 - v1, point - v1)
            c2 = glm.cross(v0 - v2, point - v2)

            if glm.dot(c0, normal) >= 0 and glm.dot(c1, normal) >= 0 and glm.dot(c2, normal) >= 0:
                return True

        return False

    def __repr__(self):
        return f"Mesh({self.name})"
