#include "gray_codes.hpp"

int main() {
    gray_codes::GrayCodeGenerator<uint8_t> generator(5, true);
    generator.generate();
}
