#include <main.hpp>

int main() {
    generateGrayCodes(5, true);
}

bool generateGrayCodes(char wordLengthBits, bool beckett) {
    Word wordLengthBytes = (wordLengthBits / 8) + (wordLengthBits % 8 != 0 ? 1 : 0);
    if (sizeof(Word) < wordLengthBytes) {
        fprintf(
            stderr,
            "Requested word length requires %d bytes, but the Word datatype is %d bytes\n",
            wordLengthBytes, sizeof(Word)
        );
        return false;
    }

    Word sequenceLength = 1 << wordLengthBits;
    std::stack<GrayCode> stack = std::stack<GrayCode>();
    std::stack<GrayCode> treeLayer = std::stack<GrayCode>();
    unsigned long count = 0;

    // If we only accept cyclic sequences, all legal codes with nonzero first words
    //  can be found by rotating one of the codes with zero as the first word,
    //  so we only need to explore the subtree starting with zero
    stack.push(GrayCode(sequenceLength, wordLengthBits));
    stack.top().available.toggle(0);
    stack.top().sequence.push_back(0);

    std::vector<std::vector<Word>> finalCodes = std::vector<std::vector<Word>>();
    std::vector<Word> curPermutation = std::vector<Word>(wordLengthBits, 0);
    for (Word i = 0; i < curPermutation.size(); i++) {
        curPermutation[i] = i;
    }
    std::vector<std::vector<Word>> bitPositionPermutations = std::vector<std::vector<Word>>();
    for (std::vector<bool>::size_type i = 0; i < factorial(wordLengthBits); i++) {
        std::next_permutation(curPermutation.begin(), curPermutation.end());
        bitPositionPermutations.push_back(std::vector(curPermutation));
    }

    while (stack.size() > 0) {
        GrayCode prevCode = stack.top();
        stack.pop();
        if (prevCode.sequence.size() == sequenceLength) {
            if ((prevCode.sequence.back() & (prevCode.sequence.back() - 1)) == 0) {
                // The last element is a single bit flip away from zero

                // Ignore codes that are isomorphic via having their words reversed + bits permuted
                bool foundReverseIsomorph = false;
                for (std::vector<Word> sequence : finalCodes) {
                    if (isIsomorphic(prevCode.sequence, sequence,
                            bitPositionPermutations, wordLengthBits)) {
                        foundReverseIsomorph = true;
                        break;
                    }
                }
                if (!foundReverseIsomorph) {
                    count++;
                    fprintf(stdout, "count = %d\n", count);
                    finalCodes.push_back(prevCode.sequence);

                    Word prevWord = 0;
                    Word oneHot = 0;
                    Word transitionIndex = 0;
                    for (Word w = 1; w < prevCode.sequence.size(); w++) {
                        oneHot = prevWord ^ prevCode.sequence[w];
                        transitionIndex = 0;
                        while ((oneHot >>= 1) > 0) transitionIndex++;
                        fprintf(stdout, "%u", transitionIndex);
                        prevWord = prevCode.sequence[w];
                    }
                    oneHot = prevWord ^ 0x00;  // get transition to 0 (xor is a nop here, silly)
                    transitionIndex = 0;
                    while ((oneHot >>= 1) > 0) transitionIndex++;
                    fprintf(stdout, "%u", transitionIndex);
                    fprintf(stdout, "\n");

                    // prevCode.print();
                }

                // if (count % 100 == 0) {
                //     fprintf(stdout, "%d\n", count);
                // }
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
                //     std::fprintf(stdout, "Denying Beckett candidate 0b");
                //     for (int i = wordLengthBits - 1; i >= 0; i--) {
                //         std::fprintf(stdout, "%s", (candidateWord & (1 << i)) > 0 ? "1" : "0");
                //     }
                //     std::fprintf(stdout, " w/times ");
                //     for (int i = wordLengthBits - 1; i >= 0; i--) {
                //         std::fprintf(stdout, "%d, ", prevCode.setTimes[i]);
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

    return true;
}

BitVector::BitVector(std::vector<Cell>::size_type numBits, bool initialValue) {
    data = std::vector<Cell>(
        (numBits / sizeof(Cell)) + (numBits % sizeof(Cell) != 0 ? 1 : 0),
        initialValue ? ~0 : 0
    );
    length = numBits;
}

void BitVector::toggle(std::vector<Cell>::size_type index) {
    data[index / sizeof(Cell)] ^= 1 << (index % sizeof(Cell));
}

bool BitVector::get(std::vector<Cell>::size_type index) {
    return (data[index / sizeof(Cell)] & (0x01 << (index % sizeof(Cell)))) > 0;
}

std::vector<BitVector::Cell>::size_type BitVector::size() {
    return length;
}

void GrayCode::openBit(Word bitIndex) {
    isOpen |= 0x01 << bitIndex;
}

bool GrayCode::isBitOpen(Word bitIndex) {
    return (isOpen & (0x01 << bitIndex)) > 0;
}

bool isIsomorphic(std::vector<Word> sequence1, std::vector<Word> sequence2,
        std::vector<std::vector<Word>> bitPositionPermutations, Word wordLengthBits) {
    bool candidateIsIsomorphic = false;
    for (std::vector<Word> permutation : bitPositionPermutations) {
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

bool getBit(Word word, Word bitIndex) {
    return (word & (0x01 << bitIndex)) > 0;
}

bool GrayCode::isOldestSetBit(Word bitIndex) {
    Word bitSetTime = setTimes[bitIndex];
    for (Word i = 0; i < setTimes.size(); i++) {
        if (i != bitIndex && (sequence.back() & (0x01 << i)) > 0 && setTimes[i] < bitSetTime) {
            return false;
        }
    }
    return true;
}

void GrayCode::print() {
    std::vector<std::string> strs = std::vector<std::string>(setTimes.size(), "");
    for (Word word : sequence) {
        for (Word i = 0; i < strs.size(); i++) {
            strs[i] += (word & 0x01) == 0x01 ? '1' : '0';
            word >>= 1;
        }
    }
    for (std::string str : strs) {
        fprintf(stdout, "%s\n", str.c_str());
    }
    fprintf(stdout, "\n");
}

std::vector<bool>::size_type factorial(Word n) {
    std::vector<bool>::size_type result = 1;
    for (Word i = n; i >= 1; i--) {
        result *= i;
    }
    return result;
}
