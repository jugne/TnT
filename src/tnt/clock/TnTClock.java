package tnt.clock;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import starbeast2.SpeciesTreeRates;
import tnt.distribution.GeneTreeIntervals;

public class TnTClock extends BranchRateModel.Base {
	public Input<GeneTreeIntervals> intervalsInput = new Input<>("geneTreeIntervals",
			"The gene tree this relaxed clock is associated with.", Input.Validate.REQUIRED);
	public Input<SpeciesTreeRates> transmissionTreeRatesInput = new Input<>("transmissionTreeRates",
			"The per-branch rates for the species tree", Input.Validate.REQUIRED);

    private int geneNodeCount;
    private double[] branchRates;
    private double[] storedBranchRates;
    private boolean needsUpdate;

    RealParameter meanRate;
    SpeciesTreeRates speciesTreeRatesX;
	GeneTreeIntervals intervals;
    
    @Override
    public void initAndValidate() {
        meanRate = meanRateInput.get();
		speciesTreeRatesX = transmissionTreeRatesInput.get();
		intervals = intervalsInput.get();
		geneNodeCount = intervals.geneTreeInput.get().getNodeCount();

        branchRates = new double[geneNodeCount];
        storedBranchRates = new double[geneNodeCount];
        needsUpdate = true;
    }

    @Override
    public boolean requiresRecalculation() {
		needsUpdate = intervalsInput.isDirty() || transmissionTreeRatesInput.isDirty() || meanRateInput.isDirty();
        return needsUpdate;
    }

    @Override
    public void store() {
        System.arraycopy(branchRates, 0, storedBranchRates, 0, branchRates.length);
        super.store();
    }

    @Override
    public void restore() {
        double[] tmpRatesArray = branchRates;
        branchRates = storedBranchRates;
        storedBranchRates = tmpRatesArray;
        super.restore();
    }

    private void update() {
        final double geneTreeRate = meanRate.getValue();
		final double[] transmissionTreeRates = speciesTreeRatesX.getRatesArray();
		final double[] trNodeOccupancy = intervals.getTrNodeOccupancy();

		final int speciesNodeCount = transmissionTreeRates.length;
        for (int i = 0; i < geneNodeCount - 1; i++) {
            double weightedSum = 0.0;
            double branchLength = 0.0;
            for (int j = 0; j < speciesNodeCount; j++) {
                // System.out.println(String.format("%d, %d: %f, %f", i, j, speciesTreeRates[j], speciesOccupancy[i * speciesNodeCount + j]));
				weightedSum += transmissionTreeRates[j] * trNodeOccupancy[i * speciesNodeCount + j];
				branchLength += trNodeOccupancy[i * speciesNodeCount + j];
            }

            branchRates[i] = geneTreeRate * weightedSum / branchLength;
			if (Double.isNaN(branchRates[i]) && branchLength == 0)
				branchRates[i] = 0;
        }
        // set the rate for the root branch of this gene to equal the input mean rate
        branchRates[geneNodeCount - 1] = geneTreeRate;

        needsUpdate = false;
    }

    @Override
    public double getRateForBranch(Node node) {
        if (needsUpdate) {
            update();
        }

        return branchRates[node.getNr()];
    }
}

