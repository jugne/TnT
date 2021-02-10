package tnt.util;

// COPIED from starbeast2 package on 18.11.2020!!

// store and return a single double value
// value if never set() is positive infinity
// if set() is called multiple times, the smallest value will be stored
public final class MinimumDouble {
    private double storedDouble;

    public MinimumDouble() {
        storedDouble = Double.POSITIVE_INFINITY;
    }

    public void set(double inputDouble) {
        if (inputDouble < storedDouble) {
            storedDouble = inputDouble;
        }
    }

    public double get() {
        return storedDouble;
    }
}
