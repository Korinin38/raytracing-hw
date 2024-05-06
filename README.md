Repository for Raytracing project.

![1024 samples per pixel](examples/sponza.png "Sponza")

## Requirements
* CMake v3.27 or newer
* GNU GCC 10.5 or newer

Optional: OpenMP

## Installation and execution
To build: [build.sh](build.sh)  
To run: [run.sh](run.sh)

To clean: [clean.sh](clean.sh)

---

Run arguments:
`run.sh input.gltf width height sample_count output.ppm`
1. `input.gltf` - Input file in [GLTF](https://registry.khronos.org/glTF/specs/2.0/glTF-2.0.html) format
2. `widht` - Width of the output image
3. `height` - Height of the output image
4. `samples` - Sample count per pixel
5. `output.ppm` - Output image in [PPM](https://ru.wikipedia.org/wiki/Portable_anymap) format
