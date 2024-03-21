import math

import numpy as np
from tqdm import tqdm
import itertools
import taichi as ti
import taichi.math as tm
import geometry as geom
import helperclasses as hc
import glm

# Ported from C++ by Melissa Katz
# Adapted from code by Loïc Nassif and Paul Kry

shadow_epsilon = 10**(-6)


@ti.data_oriented
class Scene:

    def __init__(self,
                 width: int,
                 height: int,
                 jitter: bool,
                 samples: int,
                 position: tm.vec3,
                 lookat: tm.vec3,
                 up: tm.vec3,
                 fov: float,
                 ambient: tm.vec3,
                 lights: list[hc.Light],
                 materials: list[hc.Material],
                 objects: list[geom.Geometry]
                 ):
        self.width = width  # width of image
        self.height = height  # height of image
        self.aspect = width / height  # aspect ratio
        self.jitter = jitter  # should rays be jittered
        self.samples = samples  # number of rays per pixel
        self.position = position  # camera position in 3D
        self.lookat = lookat  # camera look at vector
        self.up = up  # camera up position
        self.fov = fov  # camera field of view
        self.ambient = ambient  # ambient lighting
        self.lights = lights  # all lights in the scene
        self.materials = materials  # all materials of objects in the scene
        self.objects = ti.field(objects)  # all objects in the scene

    def render(self):
        # image = np.zeros((self.width, self.height, 3))
        image = ti.field(ti.f32, (self.width, self.height, 3))
        image.fill(0)

        cam_dir = self.position - self.lookat
        d = 1.0
        top = d * math.tan(0.5 * math.pi * self.fov / 180)
        right = self.aspect * top
        bottom = -top
        left = -right

        cam_dir_glm = glm.vec3(cam_dir.x, cam_dir.y, cam_dir.z)
        up_glm = glm.vec3(self.up.x, self.up.y, self.up.z)
        w = glm.normalize(cam_dir_glm)
        u = glm.cross(up_glm, w)
        u = glm.normalize(u)
        v = glm.cross(w, u)

        u = tm.vec3(u.x, u.y, u.z)
        v = tm.vec3(v.x, v.y, v.z)
        w = tm.vec3(w.x, w.y, w.z)

        for i, j in tqdm(
            itertools.product(range(self.width), range(self.height)),
            total=self.width * self.height,
            desc="Rendering",
        ):
            self.cast_ray(image, u, v, w, d, left, right, bottom, top, i, j)

        return self.image

    @ti.kernel
    def cast_ray(self,
                 image: ti.template(),
                 u: tm.vec3,
                 v: tm.vec3,
                 w: tm.vec3,
                 d: float,
                 left: float,
                 right: float,
                 top: float,
                 bottom: float,
                 i: int,
                 j: int
                 ):
        colour = tm.vec3(0, 0, 0)

        # TODO: Generate rays
        x = left + (right - left) * (i + 0.5) / self.width
        y = top - (top - bottom) * (j + 0.5) / self.height
        origin = self.position
        direction = tm.normalize(x * u + y * v - d * w)
        ray = hc.Ray(origin, direction)

        # TODO: Test for intersection
        intersections = []
        for obj in self.objects:
            intersect = obj.intersect(ray, hc.Intersection.default())
            if intersect is not None and intersect.time < float("inf"):
                intersections.append(intersect)

        if len(intersections) == 0:
            self.image[i, j, 0] = 0.0
            self.image[i, j, 1] = 0.0
            self.image[i, j, 2] = 0.0
            return

        intersections.sort(key=lambda x: x.time)
        firstIntersection = intersections[0]

        # colour = firstIntersection.mat.diffuse * self.ambient
        colour = tm.vec3(1, 0, 0)

        # TODO: Perform shading computations on the intersection point

        self.image[i, j, 0] = max(0.0, min(1.0, colour.x))
        self.image[i, j, 1] = max(0.0, min(1.0, colour.y))
        self.image[i, j, 2] = max(0.0, min(1.0, colour.z))
