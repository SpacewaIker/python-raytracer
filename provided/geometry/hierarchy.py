from geometry import Geometry, Intersection
import glm
import helperclasses as hc


class Hierarchy(Geometry):
    def __init__(self, name: str, gtype: str, materials: list[hc.Material], hierarchy_type: str, t: glm.vec3, r: glm.vec3, s: glm.vec3, speed):
        super().__init__(name, gtype, materials, speed)
        self.t = t
        self.M = glm.mat4(1.0)
        self.Minv = glm.mat4(1.0)
        self.hierarchy_type = hierarchy_type
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

    def intersect(self, ray: hc.Ray) -> list[Intersection]:
        m_inv_o = glm.vec3(self.Minv * glm.vec4(ray.origin, 1))
        m_inv_d = glm.vec3(self.Minv * glm.vec4(ray.direction, 0))
        m_ray = hc.Ray(m_inv_o, m_inv_d)

        intersections = []

        if self.hierarchy_type == "union":
            for child in self.children:
                intersections += child.intersect(m_ray)
        elif self.hierarchy_type == "intersection":
            for child in self.children:
                child_intersections = child.intersect(m_ray)
                for intersection in child_intersections:
                    if all([c.is_inside(intersection.position) or c == child for c in self.children]):
                        intersections.append(intersection)
        elif self.hierarchy_type == "difference":
            c1 = self.children[0].intersect(m_ray)
            c2 = self.children[1].intersect(m_ray)

            for intersection in c1:
                if not self.children[1].is_inside(intersection.position):
                    intersections.append(intersection)

            for intersection in c2:
                if self.children[0].is_inside(intersection.position):
                    intersection.mat = self.children[0].get_material(intersection.position)
                    intersection.normal = -intersection.normal
                    intersections.append(intersection)

        for intersect in intersections:
            if intersect.mat is None:
                intersect.mat = self.materials[0]
            intersect.position = (self.M * glm.vec4(intersect.position, 1)).xyz
            intersect.normal = glm.normalize(glm.transpose(self.Minv) * glm.vec4(intersect.normal, 0)).xyz

        return intersections

    def shadow_intersect(self, ray: hc.Ray, t_max: float) -> bool:
        m_inv_o = glm.vec3(self.Minv * glm.vec4(ray.origin, 1))
        m_inv_d = glm.vec3(self.Minv * glm.vec4(ray.direction, 0))
        m_ray = hc.Ray(m_inv_o, m_inv_d)

        if self.hierarchy_type == "union":
            for child in self.children:
                if child.shadow_intersect(m_ray, t_max):
                    return True
            return False
        elif self.hierarchy_type == "intersection":
            for child in self.children:
                if not child.shadow_intersect(m_ray, t_max):
                    return False
            return True
        elif self.hierarchy_type == "difference":
            c1 = self.children[0].intersect(m_ray)
            c2 = self.children[1].intersect(m_ray)

            for intersection in c1:
                if intersection.time > self.shadow_epsilon and not self.children[1].is_inside(intersection.position):
                    return True

            for intersection in c2:
                if intersection.time > self.shadow_epsilon and self.children[0].is_inside(intersection.position):
                    return True

            return False

        return False

    def is_inside(self, point: glm.vec3) -> bool:
        m_inv_p = glm.vec3(self.Minv * glm.vec4(point, 1))

        if self.hierarchy_type == "union":
            for child in self.children:
                if child.is_inside(m_inv_p):
                    return True
            return False
        elif self.hierarchy_type == "intersection":
            for child in self.children:
                if not child.is_inside(m_inv_p):
                    return False
            return True
        elif self.hierarchy_type == "difference":
            c1 = self.children[0].is_inside(m_inv_p)
            c2 = self.children[1].is_inside(m_inv_p)
            return c1 and not c2

        return False

    def get_material(self, point: glm.vec3) -> hc.Material:
        m_inv_p = glm.vec3(self.Minv * glm.vec4(point, 1))

        for child in self.children:
            if child.is_inside(m_inv_p):
                return child.get_material(m_inv_p)

        return None

    def set_scene(self, scene):
        self.scene = scene
        for child in self.children:
            child.set_scene(scene)

    def __repr__(self):
        return f"Hierarchy({self.name}, t: {self.t}, r: {self.r}, s: {self.s}, children: {len(self.children)})"
