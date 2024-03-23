import numpy
import scene_parser
import argparse
from PIL import Image

# Ported from C++ by Melissa Katz
# Adapted from code by Loïc Nassif and Paul Kry

parse = argparse.ArgumentParser()
parse.add_argument("--infile", type=str,
                   help="Name of json file that will define the scene")
parse.add_argument("--outfile", type=str, default="out.png",
                   help="Name of png that will contain the render")
args = parse.parse_args()

if __name__ == "__main__":

    full_scene = scene_parser.load_scene(args.infile)
    image = full_scene.render()
    image = numpy.rot90(image, k=1, axes=(0, 1))

    im = Image.fromarray((image * 255).astype(numpy.uint8))
    im.save(args.outfile)
    im.show()
