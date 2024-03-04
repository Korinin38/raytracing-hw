#include "core/scene.h"
#include "utils/timer.h"

#include <iostream>

static bool ProgressCallback(int progress, void *userData) {
    auto t = (timer *) userData;
    if (progress == 0)
        t->restart();
    FILE *f;
    std::cout << "\r   Render [";
    for (int i = 0; i < 10; i++) {
        std::cout << (progress / ((i + 1) * 10) ? "*" : " ");
    }
    std::cout << "] " << progress << "% (" << t->elapsed() << " s)" << std::flush;
    if (progress == 100) {
        std::cout <<"\n      " << t->elapsed() << " seconds (" << t->elapsed() * 1000 << " ms) elapsed." << std::endl;
    }
    return true;
}

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

    timer t;
    Scene scene(input);
    std::cout << "Scene loaded (" << t.elapsed() << " seconds)." << std::endl;

    t.restart();
//    scene.render(ProgressCallback);
    scene.render();

    scene.draw_into(output);
    std::cout << "Frame drawn into " << output << std::endl;
    return 0;
}
