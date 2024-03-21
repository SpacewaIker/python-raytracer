import math

import glm
import numpy as np
import itertools
import geometry as geom
import helperclasses as hc
from tqdm import tqdm

# Ported from C++ by Melissa Katz
# Adapted from code by Loïc Nassif and Paul Kry

shadow_epsilon = 10**(-6)


class Scene:

    def __init__(self,
                 width: int,
                 height: int,
                 jitter: bool,
                 samples: int,
                 position: glm.vec3,
                 lookat: glm.vec3,
                 up: glm.vec3,
                 fov: float,
                 ambient: glm.vec3,
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
        self.objects = objects  # all objects in the scene

    def render(self):

        image = np.zeros((self.width, self.height, 3))

        cam_dir = self.position - self.lookat
        d = 1.0
        top = d * math.tan(0.5 * math.pi * self.fov / 180)
        right = self.aspect * top
        bottom = -top
        left = -right

        w = glm.normalize(cam_dir)
        u = glm.cross(self.up, w)
        u = glm.normalize(u)
        v = glm.cross(w, u)

        for i, j in tqdm(
            itertools.product(range(self.width), range(self.height)),
            total=self.width * self.height,
            desc="Rendering",
        ):
            colour = glm.vec3(0, 0, 0)

            # TODO: Generate rays
            x = left + (right - left) * (i + 0.5) / self.width
            y = top - (top - bottom) * (j + 0.5) / self.height
            origin = self.position
            direction = glm.normalize(x * u + y * v - d * w)
            ray = hc.Ray(origin, direction)

            # TODO: Test for intersection
            intersections = []
            for obj in self.objects:
                intersect = obj.intersect(ray, hc.Intersection.default())
                if intersect is not None and intersect.time < float("inf"):
                    intersections.append(intersect)

            if len(intersections) == 0:
                image[i, j, 0] = 0.0
                image[i, j, 1] = 0.0
                image[i, j, 2] = 0.0
                continue

            intersections.sort(key=lambda x: x.time)
            first_intersect = intersections[0]

            colour = self.ambient * first_intersect.mat.diffuse

            # TODO: Perform shading computations on the intersection point

            image[i, j, 0] = max(0.0, min(1.0, colour.x))
            image[i, j, 1] = max(0.0, min(1.0, colour.y))
            image[i, j, 2] = max(0.0, min(1.0, colour.z))

        return image
