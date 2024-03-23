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

        progress = tqdm(total=self.width * self.height, desc="Rendering")

        for i in range(self.width):
            y = vc.bottom + 0.5 * dy

            for j in range(self.height):
                # if i == 127 and j == 96:
                #     breakpoint()

                origin = self.position
                light_dir = glm.normalize(x * vc.u + y * vc.v - vc.d * vc.w)
                ray = hc.Ray(origin, light_dir)

                self.cast_ray(ray, image, i, j)

                progress.update(1)

                y += dy

            x += dx

        return image

    def cast_ray(self, ray: hc.Ray, image, i: int, j: int):
        colour = glm.vec3(0, 0, 0)

        # get all intesections with all objects
        intersections = []
        for obj in self.objects:
            intersect = obj.intersect(ray)
            if intersect is not None and intersect.time < float("inf"):
                intersections.append(intersect)

        if len(intersections) == 0:
            image[i, j, 0] = 0.0
            image[i, j, 1] = 0.0
            image[i, j, 2] = 0.0
            return

        # get first intersection
        intersections.sort(key=lambda x: x.time)
        first_intersect = intersections[0]

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

        image[i, j, 0] = max(0.0, min(1.0, colour.x))
        image[i, j, 1] = max(0.0, min(1.0, colour.y))
        image[i, j, 2] = max(0.0, min(1.0, colour.z))
