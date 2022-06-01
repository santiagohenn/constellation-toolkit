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

    /* inputSet[]  ---> Input Array
        data[] ---> Temporary array to store current combination
        start & end ---> Staring and Ending indexes in inputSet[]
        index  ---> Current index in data[]
        sizeOfSubset ---> Size of a combination to be printed */
    private void combinationUtil(int[] inputSet, int n, int sizeOfSubset, int index, int[] data, int i) {
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
        if (i >= n)
            return;

        // current is included, put next at next location
        data[index] = inputSet[i];
        combinationUtil(inputSet, n, sizeOfSubset, index + 1, data, i + 1);

        // current is excluded, replace it with next (Note that
        // i+1 is passed, but index is not changed)
        combinationUtil(inputSet, n, sizeOfSubset, index, data, i + 1);
    }

    // The main function that prints all combinations of size r
    // in arr[] of size n. This function mainly uses combinationUtil()
    private void computeCombinations(int[] arr, int r) {
        // A temporary array to store all combination one by one
        int[] data = new int[r];

        // Print all combination using temporary array 'data[]'
        combinationUtil(arr, arr.length, r, 0, data, 0);
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
