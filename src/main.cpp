#include "geometry/scene.h"
#include <iostream>

int parse_args(int argc, char* argv[], std::string &input, std::string &output) {
    if (argc < 2 || argc > 3) {
        std::cerr << "Invalid arguments - " << argc << " (expected: 2)" << std::endl;
        return 1;
    }
    input = argv[1];
    output = "output.ppm";
    if (argc == 3)
        output = argv[2];
    return 0;
}

int main(int argc, char* argv[]) {
    std::string input;
    std::string output;
    if (parse_args(argc, argv, input, output))
        return 1;

    Scene sceneScene(input);
    std::cout << "Scene loaded." << std::endl;
    return 0;
}
