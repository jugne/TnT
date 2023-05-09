package tnt.logger;

import java.io.PrintStream;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.util.Utils;
import beast.core.CalculationNode;
import beast.core.Function;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Tree;

/**
 * @author for SA logger: Alexandra Gavryushkina
 * @author for rho samples and total samples: Ugne Stolz
 */

public class SATreeStatsLogger extends CalculationNode implements Loggable, Function {
	public Input<Tree> treeInput = new Input<Tree>("tree", "tree to report SA, rho samples and samples count for.",
			Input.Validate.REQUIRED);

	public Input<Boolean> reportSAInput = new Input<Boolean>("reportSA", "include number of sampled ancestors", true);

	public Input<Boolean> reportRhoSamplesInput = new Input<Boolean>("reportRhoSamples",
			"include number of rho samples", true);

	public Input<Boolean> reportNonSASamplesInput = new Input<Boolean>("reportNonSASample",
			"include the number of non-SA samples", true);

	public Input<Function> finalSampleOffsetInput = new Input<>("finalSampleOffset",
			"The difference in time between the final sample and the end of the BD process. " +
					"Will be set by the simulator.",
			Input.Validate.REQUIRED);

	public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
			"BDMM parameterization",
			Input.Validate.REQUIRED);


    @Override
    public void initAndValidate() {
		// nothing to do
    }

    @Override
    public void init(PrintStream out) {
		if (reportSAInput.get())
			out.print("SACount" + "\t");
		if (reportRhoSamplesInput.get())
			out.print("RhoSampleCount" + "\t");
		if (reportNonSASamplesInput.get())
			out.print("NonSASampleCount" + "\t");
    }

    @Override
    public void log(long nSample, PrintStream out) {
		final Tree tree = treeInput.get();
		Parameterization parameterization = parameterizationInput.get();
		int saCount = tree.getDirectAncestorNodeCount();
		if (reportSAInput.get())
			out.print(saCount + "\t");
		if (reportRhoSamplesInput.get()) {
			// Determine which, if any, of the leaf ages correspond exactly to
			// rho sampling times.
			Integer nRhoTip = 0;
			for (int nodeNr = 0; nodeNr < tree.getLeafNodeCount(); nodeNr++) {
				double nodeTime = parameterization.getNodeTime(tree.getNode(nodeNr),
						finalSampleOffsetInput.get().getArrayValue());
//	            double nodeTime = parameterization.getTotalProcessLength() - tree.getNode(nodeNr).getHeight();
				for (double rhoSampTime : parameterization.getRhoSamplingTimes()) {
					if (Utils.equalWithPrecision(rhoSampTime, nodeTime)) {
						nRhoTip += 1;
						break;
					}
				}
			}
			out.print(nRhoTip + "\t");
		}
		if (reportNonSASamplesInput.get())
			out.print(tree.getLeafNodeCount() - saCount + "\t");

    }

	@Override
	public void close(PrintStream out) {
		// nothing to do
    }

	@Override
	public int getDimension() {
		return 1;
    }

	@Override
	public double getArrayValue() {
		return treeInput.get().getDirectAncestorNodeCount();
    }

    @Override
	public double getArrayValue(int iDim) {
		return treeInput.get().getDirectAncestorNodeCount();
    }

}

