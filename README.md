# RayTracer

**A program capable of creating images of scenes while utilizing ray tracing algorithms.**

This is a program for a computer graphics class project that utilizes ray tracing algorithms to generate images of scenes as ASCII.ppm files. It does so by reading in an input file with various information about the scene. Information such as camera position, geometric object coordinates and attributes, light positions and color, etc. After reading an input file, it determines if enough valid information has been passed in to create a scene and generates an image file. 

## How to Start

The makefile allows for the compilation and generation of an executable simply called 'raytracer'. In order to run the raytracer, a command line argument must be provided of a txt input file that contains enough necessary information. For example, to run the raytracer for input file Sample1.txt, you would run the program with

```
./raytracer Sample1.txt
```

This would generate an ASCII.ppm file named Sample1.ppm. 

![](Sample1.ppm)

## Input Files
Example input txt files can be found in the repository. Each essential component is followed by a series of numbers. These typically represent a components 3d or 2d coordinates. For example, an eye/camera placed at coordinates (0,0,0) will have those three numbers proceding it in the input file. Further details about input file syntax can be found in the wiki.

## Capabilities

At its current state, the program is able to generate images of scenes containing spheres, triangles, and lights. Each pixel's rgb value is calculated based on camera orientation, lights, and a shape's color characteristics calculated using the [phong illumination equation](https://en.wikipedia.org/wiki/Phong_reflection_model). Textures can also be wrapped around objects as seen in Sample2.ppm. Texture images must be in ASCII.ppm format as well. Shadows and reflection are also capable.  

## Shortcomings

Currently, the program is only capable of using one texture at a time. I suspect it is due to the size of ASCII.ppm files and reading too many might be more than the program can currently handle. I suspect there is likely a workaround for this. Transparency was another feature I sought to implement, but had not been able to achieve. Order of components in the input file also matters. For example, a texture will only be applied to objects if it is declared before the object is. 
