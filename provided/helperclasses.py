import glm
from typing import Any

# Ported from C++ by Melissa Katz
# Adapted from code by Loïc Nassif and Paul Kry


class Ray:
    def __init__(self, o: glm.vec3, d: glm.vec3):
        self.origin = o
        self.direction = d

    def getDistance(self, point: glm.vec3):
        return glm.length(point - self.origin)

    def getPoint(self, t: float):
        return self.origin + self.direction * t

    def __repr__(self):
        return f"Ray(origin: {self.origin}, dir: {self.direction})"


class Material:
    def __init__(self, name: str, specular: glm.vec3, diffuse: glm.vec3, hardness: float, ID: int, mat_type: str = "diffuse", mat_tint: float = 0.0):
        self.name = name
        self.specular = specular
        self.diffuse = diffuse
        self.hardness = hardness
        self.ID = ID
        self.mat_type = mat_type
        self.refr_index = 1.0
        self.tint = mat_tint

    @staticmethod
    def default():
        name = "default"
        specular = diffuse = glm.vec3(0, 0, 0)
        hardness = ID = -1
        return Material(name, specular, diffuse, hardness, ID)

    def __repr__(self):
        return f"Material({self.name}, type: {self.mat_type}, specular: {self.specular}, diffuse: {self.diffuse}, hardness: {self.hardness}, ID: {self.ID})"


class Light:
    def __init__(self, ltype: str, name: str, colour: glm.vec3, vector: glm.vec3, power: float):
        self.type = ltype
        self.name = name
        self.colour = colour
        self.vector = vector
        self.power = power

    def __repr__(self):
        return f"Light({self.name}, type: {self.type}, colour: {self.colour}, vector: {self.vector}, power: {self.power})"


class AAInterval:
    def __init__(self, t1: float, t2: float, label: Any = None):
        self.start = min(t1, t2)
        self.end = max(t1, t2)
        self.label = label


class ViewportCamera:
    def __init__(self):
        self.focal_length = 1.0
        self.aperture = 0.0
        self.dof_samples = 1

    def set_viewport(self, width: int, height: int) -> "ViewportCamera":
        self.width = width
        self.height = height
        self.aspect = width / height
        return self

    def set_camera(self, position: glm.vec3, lookat: glm.vec3, up: glm.vec3, fov: float) -> "ViewportCamera":
        cam_dir = position - lookat
        self.position = position
        self.d = 1.0
        self.top = self.d * glm.tan(glm.radians(fov / 2))
        self.right = self.aspect * self.top
        self.bottom = -self.top
        self.left = -self.right

        self.w = glm.normalize(cam_dir)
        self.u = glm.normalize(glm.cross(up, self.w))
        self.v = glm.cross(self.w, self.u)

        return self

    def set_lens(self, focal_length: float, aperture: float, dof_samples: int) -> "ViewportCamera":
        self.focal_length = focal_length
        self.aperture = aperture
        self.dof_samples = dof_samples

        return self

    def set_motion(self, time: float, motion_samples: int, motion_final: int) -> "ViewportCamera":
        dt = time / motion_samples
        self.motion_times = [dt * i for i in range(motion_samples)]
        self.motion_times += [time] * motion_final

        return self
