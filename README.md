# Raytracer
The project includes the implementation of a raytracer with capabilities for rendering scenes, handling intersections with various shapes, computing lighting using Phong reflection, handling shadows, and recursive ray tracing for reflective surfaces.


## Project Overview
The provided framework loads scenes from a JSON configuration file and manages the output by generating a PNG image. The path to the configuration file should be passed to the program as an argument. Various test scenes are available in the subdirectory **data/task1/**. A successful execution generates a file named **\<config-name>.png** in the **output/task1/** directory.

## Usage

The easiest way to set up the project is to build a Docker image based on the provided Dockerfile and to use it within an interactive session. The raytracer processes scenes from a JSON configuration file. Upon successful execution, it generates a PNG file representing the rendered scene.

```shell
$ docker build -t cgtask1 .
```
To list the JSON files:
```shell
$ ls data/task1/*.json | xargs -n1 basename

bunny.json
bunny_spotlight.json
cone_of_destiny.json
cone_shadow.json
cone_shadow_spotlights.json
cube_basic.json
dragon.json
ray_cube_basic.json
ray_cube_oriented.json
ray_cubes.json
ray_cubes_transformed.json
ray_cube_translated.json
ray_plane_basic.json
ray_planes_cone.json
three_cones.json
tumble.json
tumble_shadow.json
x-wing.json
x-wing_spotlight.json
```
Finally, use following command to generate the output using files above, for example:
```shell
$ docker run --rm -it -v "./output:/app-data/data/task1/output/task1" cgtask1 <file name>.json
``` 



