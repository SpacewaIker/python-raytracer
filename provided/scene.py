# Thibaut Baguette
# 261001513
# COMP557/ECSE532 Assignment 4
# 2024-04-08

import glm
import numpy as np
import geometry as geom
import helperclasses as hc
from tqdm import tqdm

# Ported from C++ by Melissa Katz
# Adapted from code by Loïc Nassif and Paul Kry


class Scene:

    def __init__(self,
                 vc: hc.ViewportCamera,
                 jitter: bool,
                 samples: int,
                 ambient: glm.vec3,
                 lights: list[hc.Light],
                 materials: list[hc.Material],
                 objects: list[geom.Geometry]
                 ):
        self.vc = vc  # viewport camera
        self.jitter = jitter  # should rays be jittered
        self.samples = samples  # number of rays per pixel
        self.ambient = ambient  # ambient lighting
        self.lights = lights  # all lights in the scene
        self.materials = materials  # all materials of objects in the scene
        self.objects = objects  # all objects in the scene

    def render(self, subimage: int = 0, tasks: int = 1) -> np.ndarray:
        widths = np.array_split(np.arange(self.vc.width), tasks)
        width = widths[subimage]

        image = np.zeros((len(width), self.vc.height, 3))

        dx = (self.vc.right - self.vc.left) / self.vc.width
        x = self.vc.left + (0.5 + width[0]) * dx
        dy = (self.vc.top - self.vc.bottom) / self.vc.height

        progress = tqdm(total=len(width) * self.vc.height * self.samples * self.vc.dof_samples * len(self.vc.motion_times), desc="Rendering")

        for i in width:
            y = self.vc.bottom + 0.5 * dy

            for j in range(self.vc.height):
                colour = glm.vec3(0, 0, 0)

                base_ray_origin = self.vc.position
                base_ray_direction = glm.normalize(x * self.vc.u + y * self.vc.v - self.vc.d * self.vc.w)
                focal_point = base_ray_origin + self.vc.focal_length * base_ray_direction

                for dof_origin in self._sunflower_spread(self.vc.dof_samples, base_ray_origin, self.vc.aperture):
                    dof_direction = glm.normalize(focal_point - dof_origin)

                    for ray_origin in self._sunflower_spread(self.samples, dof_origin, 2 * (dx + dy)):
                        ray = hc.Ray(ray_origin, dof_direction)

                        if self.jitter:
                            noise = 0.1 * (dx + dy) * glm.normalize(glm.vec3(np.random.rand(), np.random.rand(), np.random.rand()))
                            ray.origin += noise

                        for time in self.vc.motion_times:
                            self.current_time = time
                            colour += self.cast_ray(ray)

                            progress.update(1)

                image[i - width[0], j] = colour / (self.samples * self.vc.dof_samples * len(self.vc.motion_times))

                y += dy

            x += dx

        return image

    def cast_ray(self, ray: hc.Ray, max_recursion: int = 10, in_shape=False) -> glm.vec3:
        if max_recursion == 0:
            return glm.vec3(0, 0, 0)

        # get all intesections with all objects
        intersections = []
        for obj in self.objects:
            intersections += obj.intersect(ray)

        if len(intersections) == 0:
            return glm.vec3(0, 0, 0)

        # get first intersection
        first_intersect = min(intersections, key=lambda x: x.time)

        # check for types of materials
        if first_intersect.mat.mat_type == "mirror":
            # reflect
            reflect_dir = glm.reflect(ray.direction, first_intersect.normal)
            reflect_pos = first_intersect.position + 0.01 * reflect_dir
            reflect_ray = hc.Ray(reflect_pos, reflect_dir)
            reflection = self.cast_ray(reflect_ray, max_recursion - 1)
            colour = self._compute_regular_lighting(ray, first_intersect)
            colour = colour * first_intersect.mat.tint + reflection * (1 - first_intersect.mat.tint)
        elif first_intersect.mat.mat_type == "refractive":
            refraction = self._compute_refraction(ray, first_intersect, in_shape, max_recursion)
            colour = self._compute_regular_lighting(ray, first_intersect)
            colour = colour * first_intersect.mat.tint + refraction * (1 - first_intersect.mat.tint)
        else:
            # else, diffuse
            colour = self._compute_regular_lighting(ray, first_intersect)

        cx = max(0.0, min(1.0, colour.x))
        cy = max(0.0, min(1.0, colour.y))
        cz = max(0.0, min(1.0, colour.z))
        return glm.vec3(cx, cy, cz)

    phi = (1 + np.sqrt(5)) / 2

    def _sunflower_spread(self, num_points: int, origin: glm.vec3, radius: float) -> list[glm.vec3]:
        """
        Adapted from:
        https://stackoverflow.com/questions/28567166/uniformly-distribute-x-points-inside-a-circle
        """
        points = []
        angle_stride = 2 * np.pi / self.phi

        for k in range(1, num_points + 1):
            if k > num_points:
                r = radius
            else:
                r = radius * np.sqrt(k - 0.5) / np.sqrt(num_points - 0.5)
            theta = k * angle_stride
            x = r * np.cos(theta) + origin.x
            y = r * np.sin(theta) + origin.y
            points.append(glm.vec3(x, y, origin.z))

        return points

    def _compute_regular_lighting(self, ray: hc.Ray, intersection: geom.Intersection) -> glm.vec3:
        colour = glm.vec3(0, 0, 0)

        if isinstance(intersection.geometry, geom.simple_geometry.AABB) or isinstance(intersection.geometry, geom.simple_geometry.Plane):
            diffuse = intersection.geometry.get_diffuse(intersection.position)
        else:
            diffuse = intersection.mat.diffuse

        for light in self.lights:
            # check shadow ray
            skip = False

            # create shadow ray
            if light.type == "point":
                shadow_ray = hc.Ray(intersection.position, light.vector - intersection.position)
                t_max = 1.0
            elif light.type == "directional":
                shadow_ray = hc.Ray(intersection.position, -light.vector)
                t_max = float("inf")

            # check if shadow ray intersects with any object
            for obj in self.objects:
                if obj.shadow_intersect(shadow_ray, t_max):
                    skip = True
                    break

            if skip:
                continue

            # compute light direction
            if light.type == "point":
                light_dir = glm.normalize(light.vector - intersection.position)
            elif light.type == "directional":
                light_dir = glm.normalize(-light.vector)

            # compute lighting
            normal = intersection.normal

            lambert = diffuse * max(0.0, glm.dot(normal, light_dir))

            half_vect = glm.normalize(light_dir - ray.direction)
            specular = intersection.mat.specular * max(0.0, glm.dot(normal, half_vect)) ** intersection.mat.hardness

            colour += light.colour * light.power * (lambert + specular)

        colour += self.ambient * diffuse

        return colour

    def _compute_refraction(self, ray: hc.Ray, intersection: geom.Intersection, in_shape, max_recursion) -> glm.vec3:
        if in_shape:
            eta = intersection.mat.refr_index
            intersection.normal = -intersection.normal
        else:
            eta = 1.0 / intersection.mat.refr_index

        refract_dir = glm.refract(ray.direction, intersection.normal, eta)

        if refract_dir == glm.vec3(0, 0, 0):  # total internal reflection
            return glm.vec3(0, 0, 0)
            reflect_dir = glm.reflect(ray.direction, intersection.normal)
            reflect_pos = intersection.position + 0.001 * reflect_dir
            reflect_ray = hc.Ray(reflect_pos, reflect_dir)
            reflection = self.cast_ray(reflect_ray, max_recursion - 1, in_shape)
            return reflection

        refract_pos = intersection.position + 0.0001 * refract_dir
        refract_ray = hc.Ray(refract_pos, refract_dir)
        refraction = self.cast_ray(refract_ray, max_recursion - 1, not in_shape)
        return refraction
