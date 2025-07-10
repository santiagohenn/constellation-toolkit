package constellation.tools.math;

import java.util.ArrayList;
import java.util.List;

public class Combination {

    public List<List<Integer>> allCombinations = new ArrayList<>();
    public int numberOfMembers;
    public int maximumSizeOfSubset;

    public Combination() {

    }

    public Combination(int numberOfMembers, int maximumSizeOfSubset) {
        this.numberOfMembers = numberOfMembers;
        this.maximumSizeOfSubset = maximumSizeOfSubset;
    }

    /**
     * Generates all combinations of a given size from the input set.
     *
     * @param inputSet          the input array
     * @param inputArraySize    the size of the input array
     * @param sizeOfSubset  the size of each combination to be generated
     * @param index         the current index in the temporary array
     * @param data          the temporary array to store the current combination
     * @param i             the current index in the input array
     */
    private void combinationUtil(int[] inputSet, int inputArraySize, int sizeOfSubset, int index, int[] data, int i) {
        // Current combination is ready to be saved, save it
        if (index == sizeOfSubset) {
            List<Integer> subset = new ArrayList<>();

            for (Integer member : data) {
                subset.add(member);
            }

            allCombinations.add(subset);

            return;
        }

        // When no more elements are there to put in data[]
        if (i >= inputArraySize)
            return;

        // current is included, put next at next location
        data[index] = inputSet[i];
        combinationUtil(inputSet, inputArraySize, sizeOfSubset, index + 1, data, i + 1);

        // current is excluded, replace it with next (Note that
        // i+1 is passed, but index is not changed)
        combinationUtil(inputSet, inputArraySize, sizeOfSubset, index, data, i + 1);
    }

    /**
     * The main function that computes all combinations of a specified size from the given array.
     * This function uses combinationUtil() to generate the combinations.
     *
     * @param inputArray the input array from which combinations are to be generated
     * @param sizeOfEachCombination   the size of each combination
     */
    private void computeCombinations(int[] inputArray, int sizeOfEachCombination) {
        int[] data = new int[sizeOfEachCombination];
        combinationUtil(inputArray, inputArray.length, sizeOfEachCombination, 0, data, 0);
    }

    public List<List<Integer>> computeCombinations() {
        return computeCombinations(this.numberOfMembers);
    }

    public List<List<Integer>> computeCombinations(int numberOfMembers) {

        int[] inputSet = new int[numberOfMembers];

        for (int member = 0; member < inputSet.length; member++) {
            inputSet[member] = member;
        }

        for (int sizeOfSubset = 1; sizeOfSubset <= maximumSizeOfSubset; sizeOfSubset++) {
            computeCombinations(inputSet, sizeOfSubset);
        }

        return this.allCombinations;

    }

    public List<List<Integer>> computeCombinations(List<Integer> allMembers, int sizeOfSubset) {

        allCombinations.clear();

        // A temporary array to store all combination one by one
        int[] data = new int[sizeOfSubset];
        int[] inputSet = new int[allMembers.size()];

        int index = 0;
        for (Integer member : allMembers) {
            inputSet[index++] = member;
        }

        combinationUtil(inputSet, allMembers.size(), sizeOfSubset, 0, data, 0);

        return this.allCombinations;

    }
}
