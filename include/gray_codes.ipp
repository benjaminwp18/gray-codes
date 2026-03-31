#pragma once

#include <cstddef>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <algorithm>
#include <stdexcept>

namespace gray_codes {

template <typename Word>
GrayCodeGenerator<Word>::GrayCodeGenerator(Word wordLengthBits, bool beckett) :
        wordLengthBits(wordLengthBits), beckett(beckett) {

    Word wordLengthBytes = (wordLengthBits / 8) + (wordLengthBits % 8 != 0 ? 1 : 0);
    if (sizeof(Word) < wordLengthBytes) {
        throw std::invalid_argument(
            "Requested word length requires " + std::to_string(wordLengthBytes) +
            " bytes, but the Word datatype is " + std::to_string(sizeof(Word)) + " bytes\n"
        );
    }

    sequenceLength = 1 << wordLengthBits;

    // TODO: optionally use bit opening lexographic guarantees to generate
    //  correct permutations on the fly instead of looping through
    //  bitPositionPermutations for every evaluation
    std::vector<Word> curPermutation = std::vector<Word>(wordLengthBits, 0);
    for (Word i = 0; i < curPermutation.size(); i++) {
        curPermutation[i] = i;
    }
    bitPositionPermutations = std::vector<std::vector<Word>>();
    for (std::vector<bool>::size_type i = 0; i < factorial(wordLengthBits); i++) {
        std::next_permutation(curPermutation.begin(), curPermutation.end());
        bitPositionPermutations.push_back(std::vector(curPermutation));
    }
}

template <typename Word>
std::vector<typename GrayCodeGenerator<Word>::GrayCode>
GrayCodeGenerator<Word>::generate(GrayCode stub, Word maxDepth) {
    std::stack<GrayCode> stack = std::stack<GrayCode>();
    std::stack<GrayCode> treeLayer = std::stack<GrayCode>();
    std::vector<GrayCode> finalCodes = std::vector<GrayCode>();
    unsigned long count = 0;

    if (stub.sequence.size() > 0) {
        stack.push(stub);
    }
    else {
        // If we only accept cyclic sequences, all legal codes with nonzero first words
        //  can be found by rotating one of the codes with zero as the first word,
        //  so we only need to explore the subtree starting with zero
        stack.push(GrayCode(sequenceLength, wordLengthBits));
        stack.top().available.toggle(0);
        stack.top().sequence.push_back(0);
    }

    while (stack.size() > 0) {
        GrayCode prevCode = stack.top();
        stack.pop();
        if (maxDepth > 0 && prevCode.sequence.size() == maxDepth) {
            finalCodes.push_back(prevCode);
        }
        else if (prevCode.sequence.size() == sequenceLength) {
            if ((prevCode.sequence.back() & (prevCode.sequence.back() - 1)) == 0) {
                // The last element is a single bit flip away from zero
                insertIfUnique(prevCode, finalCodes);
            }
        }
        else {
            bool openedBit = false;
            for (Word b = 0; b < wordLengthBits; b++) {
                Word candidateWord = prevCode.sequence.back() ^ (0x01 << b);
                if (prevCode.available.get(candidateWord) &&
                    // For Beckett-Gray codes: setting a bit or unsetting the earliest set bit
                    (!beckett || candidateWord > prevCode.sequence.back() ||
                     prevCode.isOldestSetBit(b)) &&
                    // The bit has already been modified in this sequence or
                    //  we haven't opened a bit yet at this tree layer
                    (!openedBit || prevCode.isBitOpen(b)) &&
                    // If there's a candidate at a pendant vertex,
                    //  make sure this word is a pendant
                    (!prevCode.nextWordIsPendant ||
                     prevCode.numAvailableNeighbors[candidateWord] == 1)) {

                    GrayCode newCode = GrayCode(prevCode);
                    newCode.available.toggle(candidateWord);
                    newCode.sequence.push_back(candidateWord);
                    newCode.nextWordIsPendant = false;
                    if (candidateWord > prevCode.sequence.back()) {
                        newCode.setTimes[b] = prevCode.sequence.size();
                    }

                    if (!newCode.isBitOpen(b)) {
                        // "Open" the bit if this is the first time it's been modified.
                        // It's now unique. Modifying it in the same tree layer as other
                        //  bits won't create bit-rearrange isomorphisms.
                        newCode.openBit(b);
                        openedBit = true;
                    }

                    // TODO: check for bridges between Hamming classes ala Sawada & Wong
                    // Consider counting undirected edges or nodes in each class
                    // Consider using the pendant neighbor loop to get an accurate
                    //  count of the edges connected to this node & update edge
                    //  counts more accurately (O(n) instead of O(1))

                    // Hamilton cycle pendant vertices heuristic (Sawada & Wong)
                    Word unvisitedPendantNeighbors = 0;
                    Word neighbor;
                    for (Word deltaBit = 0; deltaBit < wordLengthBits; deltaBit++) {
                        neighbor = candidateWord ^ (0x01 << deltaBit);
                        if (newCode.numAvailableNeighbors[neighbor]-- <= 1
                                && newCode.available.get(neighbor)) {
                            unvisitedPendantNeighbors++;
                            // If neighbor is a pendant vertex, we have to visit it next
                            newCode.nextWordIsPendant = true;
                            if (unvisitedPendantNeighbors >= 2) {
                                break;
                            }
                        }
                    }
                    if (unvisitedPendantNeighbors >= 2) {
                        // If there exist >=2 pendant vertices that we haven't visited yet
                        //  there's no way to complete the Hamilton cycle, so this prefix
                        //  can't be part of a Gray code
                        continue;
                    }

                    treeLayer.push(newCode);

                    if (prevCode.nextWordIsPendant) {
                        // If this was a forced pendant vertex addition, we can't try
                        //  any more words
                        break;
                    }
                }
                // else {
                //     printf("Denying Beckett candidate 0b");
                //     for (int i = wordLengthBits - 1; i >= 0; i--) {
                //         printf("%s", (candidateWord & (1 << i)) > 0 ? "1" : "0");
                //     }
                //     printf(" w/times ");
                //     for (int i = wordLengthBits - 1; i >= 0; i--) {
                //         printf("%d, ", prevCode.setTimes[i]);
                //     }
                //     std::fprintf(
                //         stdout,
                //         " setting a bit: %s, oldest bit: %s, bit open: %s\n",
                //         candidateWord > prevCode.sequence.back() ? "true" : "false",
                //         prevCode.isOldestSetBit(b) ? "true" : "false",
                //         prevCode.isBitOpen(b) ? "true" : "false"
                //     );
                //     prevCode.print();
                //     // prevCode.debugIsOldestSetBit(b);
                // }
            }
            transferStack(treeLayer, stack);
        }
    }

    return finalCodes;
}

template <typename Word>
bool GrayCodeGenerator<Word>::isIsomorphic(std::vector<Word> sequence1,
                                           std::vector<Word> sequence2) {
    bool candidateIsIsomorphic = false;
    for (std::vector<Word> &permutation : bitPositionPermutations) {
        candidateIsIsomorphic = true;
        // Cyclic codes are reversed and rotated by 1 because they can't end in 0
        for (Word w = 1; w < sequence1.size(); w++) {
            for (Word b = 0; b < wordLengthBits; b++) {
                if (getBit(sequence1[w], b) !=
                    getBit(sequence2[sequence2.size() - w], permutation[b])) {
                    // std::fprintf(
                    //     stdout, "Comparing bit %u, %u (%u) to %u, %u (%u)\n",
                    //     w, b,
                    //     getBit(sequence1[w], b),
                    //     sequence2.size() - w, permutation[b],
                    //     getBit(sequence2[sequence2.size() - w], permutation[b])
                    // );
                    candidateIsIsomorphic = false;
                    break;
                }
            }
            if (!candidateIsIsomorphic) break;
        }
        if (candidateIsIsomorphic) {
            return true;
        }
    }
    return false;
}

template <typename Word>
void GrayCodeGenerator<Word>::insertIfUnique(GrayCode newCode, std::vector<GrayCode> &codes) {
    // Ignore codes that are isomorphic via having their words reversed + bits permuted
    bool foundReverseIsomorph = false;
    for (GrayCode &savedCode : codes) {
        if (isIsomorphic(newCode.sequence, savedCode.sequence)) {
            foundReverseIsomorph = true;
            break;
        }
    }
    if (!foundReverseIsomorph) {
        codes.push_back(newCode);
    }
}


/* ========== GrayCode ============= */

template <typename Word>
GrayCodeGenerator<Word>::GrayCode::GrayCode(Word sequenceLength,
                                            Word wordLengthBits) :
    available(sequenceLength, true), sequence(),
    setTimes(wordLengthBits, 0), isOpen(0),
    numAvailableNeighbors(sequenceLength, wordLengthBits),
    nextWordIsPendant(false) {
    sequence.reserve(sequenceLength);
}

template <typename Word>
void GrayCodeGenerator<Word>::GrayCode::openBit(Word bitIndex) {
    isOpen |= 0x01 << bitIndex;
}
template <typename Word>
bool GrayCodeGenerator<Word>::GrayCode::isBitOpen(Word bitIndex) {
    return (isOpen & (0x01 << bitIndex)) > 0;
}

template <typename Word>
bool GrayCodeGenerator<Word>::GrayCode::isOldestSetBit(Word bitIndex) {
    Word bitSetTime = setTimes[bitIndex];
    // TODO: use wordLengthBits instead of setTimes.size() etc.
    for (Word i = 0; i < setTimes.size(); i++) {
        if (i != bitIndex && (sequence.back() & (0x01 << i)) > 0 && setTimes[i] < bitSetTime) {
            return false;
        }
    }
    return true;
}

template <typename Word>
void GrayCodeGenerator<Word>::GrayCode::print() {
    std::vector<std::string> strs = std::vector<std::string>(setTimes.size(), "");
    for (Word &word : sequence) {
        for (Word i = 0; i < strs.size(); i++) {
            strs[i] += (word & 0x01) == 0x01 ? '1' : '0';
            word >>= 1;
        }
    }
    for (std::string &str : strs) {
        printf("%s\n", str.c_str());
    }
    printf("\n");
}


/* ============= BitVector ============= */

BitVector::BitVector(VectorSizeType numBits, bool initialValue) {
    data = std::vector<Cell>(
        (numBits / sizeof(Cell)) + (numBits % sizeof(Cell) != 0 ? 1 : 0),
        initialValue ? ~0 : 0
    );
    length = numBits;
}

void BitVector::toggle(VectorSizeType index) {
    data[index / sizeof(Cell)] ^= 1 << (index % sizeof(Cell));
}

bool BitVector::get(VectorSizeType index) {
    return (data[index / sizeof(Cell)] & (0x01 << (index % sizeof(Cell)))) > 0;
}

VectorSizeType BitVector::size() {
    return length;
}


/* ============= STATIC METHODS ============= */

template <typename Word>
bool GrayCodeGenerator<Word>::getBit(Word word, Word bitIndex) {
    return (word & (0x01 << bitIndex)) > 0;
}

template <typename Word>
std::vector<bool>::size_type GrayCodeGenerator<Word>::factorial(Word n) {
    std::vector<bool>::size_type result = 1;
    for (Word i = n; i >= 1; i--) {
        result *= i;
    }
    return result;
}

namespace {
    template <typename T>
    void transferStack(std::stack<T> &srcStack, std::stack<T> &dstStack) {
        while (!srcStack.empty()) {
            dstStack.push(srcStack.top());
            srcStack.pop();
        }
    }
}

}  // namespace gray_codes
