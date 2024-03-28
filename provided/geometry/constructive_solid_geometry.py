from geometry import Geometry, Intersection
import helperclasses as hc
import glm


class CSGUnion(Geometry):
    def __init__(self, name: str, materials: list[hc.Material], children: list[Geometry]):
        super().__init__(name, "union", materials)
        self.children = children

    def intersect(self, ray: hc.Ray) -> list[Intersection]:
        intersections = []

        for child in self.children:
            intersections += child.intersect(ray)

        return intersections

    def shadow_intersect(self, ray: hc.Ray, t_max: float) -> bool:
        for child in self.children:
            if child.shadow_intersect(ray, t_max):
                return True

        return False

    def is_inside(self, point: glm.vec3) -> bool:
        for child in self.children:
            if child.is_inside(point):
                return True

        return False

    def get_material(self, point: glm.vec3) -> hc.Material:
        for child in self.children:
            if child.is_inside(point):
                return child.get_material(point)

        return None

    def __repr__(self):
        return f"CSGUnion({self.children})"


class CSGIntersection(Geometry):
    def __init__(self, name: str, materials: list[hc.Material], children: list[Geometry]):
        super().__init__(name, "intersection", materials)
        self.children = children

    def intersect(self, ray: hc.Ray) -> list[Intersection]:
        intersections = []

        for child in self.children:
            child_intersections = child.intersect(ray)
            for intersection in child_intersections:
                if all([c.is_inside(intersection.position) or c == child for c in self.children]):
                    intersections.append(intersection)

        return intersections

    def shadow_intersect(self, ray: hc.Ray, t_max: float) -> bool:
        for child in self.children:
            if not child.shadow_intersect(ray, t_max):
                return False
        return True

    def is_inside(self, point: glm.vec3) -> bool:
        for child in self.children:
            if not child.is_inside(point):
                return False
        return True

    def get_material(self, point: glm.vec3) -> hc.Material:
        for child in self.children:
            if child.is_inside(point):
                return child.get_material(point)

        return None

    def __repr__(self):
        return f"CSGIntersection({self.children})"


class CSGDifference(Geometry):
    def __init__(self, name: str, materials: list[hc.Material], child1: Geometry, child2: Geometry):
        super().__init__(name, "difference", materials)
        self.child1 = child1
        self.child2 = child2

    def intersect(self, ray: hc.Ray) -> list[Intersection]:
        c1 = self.child1.intersect(ray)
        c2 = self.child2.intersect(ray)

        intersections = []

        for intersection in c1:
            if not self.child2.is_inside(intersection.position):
                intersections.append(intersection)

        for intersection in c2:
            if self.child1.is_inside(intersection.position):
                intersection.mat = self.child1.get_material(intersection.position)
                intersection.normal = -intersection.normal
                intersections.append(intersection)

        return intersections

    def shadow_intersect(self, ray: hc.Ray, t_max: float) -> bool:
        c1 = self.child1.intersect(ray)
        c2 = self.child2.intersect(ray)

        intersections = list(filter(lambda i: i.time > self.shadow_epsilon, c1 + c2))
        intersections.sort(key=lambda x: x.time)

        if len(intersections) == 0:
            return False

        return intersections[0] in c1 or intersections[-1] in c1

    def is_inside(self, point: glm.vec3) -> bool:
        c1 = self.child1.is_inside(point)
        c2 = self.child2.is_inside(point)
        return c1 and not c2

    def get_material(self, point: glm.vec3) -> hc.Material:
        if self.child1.is_inside(point):
            return self.child1.get_material(point)
        return None

    def __repr__(self):
        return f"CSGDifference({self.child1}, {self.child2})"
