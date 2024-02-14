#include "core/scene.h"
#include <iostream>

void parse_args(int argc, char *argv[], std::string &input, std::string &output) {
    if (argc < 2 || argc > 3) {
        throw std::runtime_error("Invalid arguments - " + std::to_string(argc) + " (expected: 2)");
    }

    input = argv[1];
    output = "output.ppm";
    if (argc == 3)
        output = argv[2];
}

int main(int argc, char *argv[]) {
    std::string input;
    std::string output;
    parse_args(argc, argv, input, output);

    Scene scene(input);
    std::cout << "Scene loaded." << std::endl;

    scene.render();
    std::cout << "Frame rendered." << std::endl;

    scene.draw_into(output);
    return 0;
}
