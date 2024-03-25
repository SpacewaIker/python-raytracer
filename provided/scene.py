import glm
import numpy as np
import geometry as geom
import helperclasses as hc
from tqdm import tqdm

# Ported from C++ by Melissa Katz
# Adapted from code by Loïc Nassif and Paul Kry


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
        vc = hc.ViewportCamera(self.position, self.lookat, self.up, self.fov, self.aspect)

        image = np.zeros((self.width, self.height, 3))

        dx = (vc.right - vc.left) / self.width
        x = vc.left + 0.5 * dx
        dy = (vc.top - vc.bottom) / self.height

        progress = tqdm(total=self.width * self.height * self.samples, desc="Rendering")

        for i in range(self.width):
            y = vc.bottom + 0.5 * dy

            for j in range(self.height):
                colour = glm.vec3(0, 0, 0)

                for ray_origin in self._sunflower_spread(self.samples, self.position, 2 * (dx + dy)):
                    light_dir = glm.normalize(x * vc.u + y * vc.v - vc.d * vc.w)
                    ray = hc.Ray(ray_origin, light_dir)

                    if self.jitter:
                        noise = 10 ** (-3) * glm.normalize(glm.vec3(np.random.rand(), np.random.rand(), np.random.rand()))
                        ray.direction = glm.normalize(ray.direction + noise)

                    colour += self.cast_ray(ray)

                    progress.update(1)

                image[i, j] = colour / self.samples

                y += dy

            x += dx

        return image

    def cast_ray(self, ray: hc.Ray, max_recursion: int = 10) -> glm.vec3:
        if max_recursion == 0:
            return glm.vec3(0, 0, 0)

        # get all intesections with all objects
        intersections = []
        for obj in self.objects:
            intersect = obj.intersect(ray)
            if intersect is not None and intersect.time < float("inf"):
                intersections.append(intersect)

        if len(intersections) == 0:
            return glm.vec3(0, 0, 0)

        # get first intersection
        first_intersect = min(intersections, key=lambda x: x.time)

        # check for types of materials
        if first_intersect.mat.mat_type == "mirror":
            # reflect
            reflect_dir = glm.reflect(ray.direction, first_intersect.normal)
            reflect_pos = first_intersect.position + 0.001 * reflect_dir
            reflect_ray = hc.Ray(reflect_pos, reflect_dir)
            reflection = self.cast_ray(reflect_ray, max_recursion - 1)
            return reflection

        # else, diffuse

        # lighting
        colour = self.ambient * first_intersect.mat.diffuse
        for light in self.lights:
            # check shadow ray
            skip = False

            # create shadow ray
            if light.type == "point":
                shadow_ray = hc.Ray(first_intersect.position, light.vector - first_intersect.position)
                t_max = 1.0
            elif light.type == "directional":
                shadow_ray = hc.Ray(first_intersect.position, -light.vector)
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
                light_dir = glm.normalize(light.vector - first_intersect.position)
            elif light.type == "directional":
                light_dir = glm.normalize(-light.vector)

            # compute lighting
            normal = first_intersect.normal

            lambert = first_intersect.mat.diffuse * light.power * max(0.0, glm.dot(normal, light_dir))

            half_vect = glm.normalize(light_dir - ray.direction)
            specular = first_intersect.mat.specular * light.power * max(0.0, glm.dot(normal, half_vect)) ** first_intersect.mat.hardness

            colour += light.colour * (lambert + specular)

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
