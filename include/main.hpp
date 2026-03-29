#include <cstddef>
#include <vector>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <stack>
#include <cstdint>
#include <algorithm>

#ifndef _MAIN_INCLUDE
#define _MAIN_INCLUDE

// Max word size is max index size for vectors, as the sequence will be represented by
//  a vector, and a sequence must exhaust all values of the word
typedef uint8_t Word;
static_assert(sizeof(std::vector<bool>::size_type) >= sizeof(Word));

// Moves all the elements in srcStack to dstStack, reversing their order.
template <typename T> void transferStack(std::stack<T> &srcStack, std::stack<T> &dstStack) {
    while (!srcStack.empty()) {
        dstStack.push(srcStack.top());
        srcStack.pop();
    }
}
std::vector<bool>::size_type factorial(Word n);
bool getBit(Word word, Word bitIndex);
bool generateGrayCodes(char wordLengthBits, bool cyclic, bool beckett);
bool isIsomorphic(std::vector<Word> sequence1, std::vector<Word> sequence2,
    std::vector<std::vector<Word>> bitPositionPermutations, bool cyclic, Word wordLengthBits);

class BitVector {
    private:
        typedef uint8_t Cell;
        std::vector<Cell> data;
        std::vector<Cell>::size_type length;
    public:
        BitVector(std::vector<Cell>::size_type numBits, bool initialValue);
        // BitVector(const BitVector &other) : data(other.data), length(other.length) {}
        BitVector(const BitVector &other) = default;

        void toggle(std::vector<Cell>::size_type index);
        bool get(std::vector<Cell>::size_type index);
        std::vector<Cell>::size_type size();
};

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

        GrayCode(Word sequenceLength, Word wordLengthBits) :
                available(sequenceLength, true), sequence(),
                setTimes(wordLengthBits, 0), isOpen(0),
                numAvailableNeighbors(sequenceLength, wordLengthBits),
                nextWordIsPendant(false) {
            sequence.reserve(sequenceLength);
        }
        // GrayCode(const GrayCode &other) :
        //     available(other.available), sequence(other.sequence),
        //     setTimes(other.setTimes), isOpen(other.isOpen),
        //     numAvailableNeighbors(other.numAvailableNeighbors),
        //     nextWordIsPendant(other.nextWordIsPendant) {}
        GrayCode(const GrayCode &other) = default;

        // void resetOpenBits();
        bool isBitOpen(Word bitIndex);
        void openBit(Word bitIndex);

        bool isOldestSetBit(Word bitIndex);
        void print();
};

#endif
