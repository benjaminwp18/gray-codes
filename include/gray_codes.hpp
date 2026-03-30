#pragma once

#include <vector>
#include <cstdint>
#include <stack>

namespace gray_codes {

namespace {
    typedef std::vector<bool>::size_type VectorSizeType;
}

class BitVector {
    private:
        typedef uint8_t Cell;
        std::vector<Cell> data;
        VectorSizeType length;
    public:
        BitVector(VectorSizeType numBits, bool initialValue);
        BitVector(const BitVector &other) = default;

        void toggle(VectorSizeType index);
        bool get(VectorSizeType index);
        VectorSizeType size();
};

/**
 * The GrayCodeGenerator is parameterized by the binary word (string) type,
 * which must be large enough to store wordLengthBits bits and no larger than
 * the maximum vector size type, defined above as std::vector<bool>::size_type.
 * Gray codes are typically parameterized by n = wordLengthBits in the
 * literature, giving a total sequence length of 2^n words.
 */
template <typename Word>
class GrayCodeGenerator {
    static_assert(sizeof(VectorSizeType) >= sizeof(Word));

    public:
        GrayCodeGenerator(Word wordLengthBits, bool beckett);

        void generate();

        // TODO: consider making GrayCode a struct
        class GrayCode {
            public:
                // The sequence of binary words representing this Gray code
                std::vector<Word> sequence;
                // available[w] = true <=> binary word w is not this code yet
                BitVector available;
                // setTimes[b] = the index in sequence at which bit b was last set to 1
                std::vector<Word> setTimes;
                // Bit b of isOpen = 1 <=> bit b has been written to previously in this code
                Word isOpen;
                // numAvailableNeighbors[w] = # of words adjacent to w that are available
                std::vector<Word> numAvailableNeighbors;
                // True <=> one of this word's neighbors is an available pendant vertex
                bool nextWordIsPendant;

                GrayCode(Word sequenceLength, Word wordLengthBits);
                GrayCode(const GrayCode &other) = default;

                // void resetOpenBits();
                bool isBitOpen(Word bitIndex);
                void openBit(Word bitIndex);

                bool isOldestSetBit(Word bitIndex);
                void print();
        };

    private:
        Word wordLengthBits;
        Word sequenceLength;
        bool beckett;
        std::vector<std::vector<Word>> bitPositionPermutations;

        bool isIsomorphic(std::vector<Word> sequence1,
                          std::vector<Word> sequence2);

        static bool getBit(Word word, Word bitIndex);
        static VectorSizeType factorial(Word n);
};

namespace {
    template <typename T>
    void transferStack(std::stack<T> &srcStack, std::stack<T> &dstStack);
}

}  // namespace gray_codes

#include "gray_codes.ipp"
