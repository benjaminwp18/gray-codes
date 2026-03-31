#include "gray_codes.hpp"
#include <omp.h>
#include <fstream>

typedef uint8_t Word;
typedef gray_codes::GrayCodeGenerator<Word> Generator;

int main() {
    Word n = 6;
    Word maxDepth = 20;

    Generator generator(n, true);
    std::vector<Generator::GrayCode> results = std::vector<Generator::GrayCode>();
    Word completedStubs = 0;

    printf("Generating stubs to a depth of %u with n = %u\n", maxDepth, n);
    std::vector<Generator::GrayCode> stubs = generator.generateStubs(maxDepth);
    std::vector<std::vector<Generator::GrayCode>> thread_results(stubs.size());
    printf("Generated %lu stubs\n", stubs.size());

    printf("Distributing stubs among max %lu OpenMP parallel threads\n", omp_get_max_threads());

    #pragma omp parallel for default(none) shared(thread_results, generator, stubs, completedStubs, stdout)
    for (size_t i = 0; i < stubs.size(); ++i) {
        // printf("Stub %lu assigned to thread %d\n", i, omp_get_thread_num());
        thread_results[i] = generator.generate(stubs[i]);
        // printf("Thread %d finished stub %lu with %lu candidate codes\n",
        //     omp_get_thread_num(), i, thread_results[i].size());

        #pragma omp atomic
        completedStubs++;

        // if (omp_get_thread_num() == 0) {
        printf("Progress: %u / %lu (%d%%)\n", completedStubs, stubs.size(),
            (int)((float)completedStubs / stubs.size() * 100));
        fflush(stdout);
        // }
    }

    for (const std::vector<Generator::GrayCode> &local_results : thread_results) {
        // results.insert(results.end(), local_results.begin(), local_results.end());
        for (const Generator::GrayCode &code : local_results) {
            generator.insertIfUnique(code, results);
        }
    }

    FILE *resultsFile = fopen("results.log", "w");

    printf("\nFinal codes:\n");

    for (const Generator::GrayCode &result : results) {
        Word prevWord = 0;
        Word oneHot = 0;
        Word transitionIndex = 0;
        // printf("\t");
        for (Word w = 1; w < result.sequence.size(); w++) {
            oneHot = prevWord ^ result.sequence[w];
            transitionIndex = 0;
            while ((oneHot >>= 1) > 0) transitionIndex++;
            fprintf(resultsFile, "%u", transitionIndex);
            // printf("%u", transitionIndex);
            prevWord = result.sequence[w];
        }
        oneHot = prevWord ^ 0x00;  // get transition to 0 (xor is a nop here, silly)
        transitionIndex = 0;
        while ((oneHot >>= 1) > 0) transitionIndex++;
        fprintf(resultsFile, "%u\n", transitionIndex);
        // printf("%u\n", transitionIndex);
    }
    fprintf(resultsFile, "Total count = %lu for n = %u\n", results.size(), n);
    printf("Total count = %lu for n = %u\n", results.size(), n);

    fclose(resultsFile);
}
