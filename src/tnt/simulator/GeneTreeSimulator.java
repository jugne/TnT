package tnt.simulator;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import beast.core.Input;
import beast.core.Runnable;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeTraceAnalysis;
import starbeast2.utils.SimulatedGeneTree;

/**
 * @author Ugne Stolz Adapted from Tim Vaughan's GeneTreeSimulator for gene
 *         trees conditioned on species tree in StarBeast2:
 *         https://doi.org/10.1093/molbev/msx126
 */

public class GeneTreeSimulator extends Runnable {

	public Input<Tree> transmissionTreeInput = new Input<>("transmissionTree",
			"Transmission tree on which to simulate gene trees", Input.Validate.REQUIRED);

	public Input<TraitSet> sampleCountsInput = new Input<>("sampleCounts",
			"TraitSet defining number of  samples per node in species tree.", Input.Validate.REQUIRED);

	// TODO transmission history traitset

	public Input<Integer> nSimsInput = new Input<>("nSims",
			"Number of gene trees to simulate from given sample distribution.", Input.Validate.REQUIRED);

	public Input<String> fileNameInput = new Input<>("fileName", "Name of file to which gene trees will be written.",
			Input.Validate.REQUIRED);

	public Input<String> reportFileNameInput = new Input<>("reportFileName",
			"Name of file to which topology distribution report will be written.");

	public Input<Double> credibilityThresholdInput = new Input<>("credibilityThreshold",
			"Maximum probability of topologies included in credible set written to report file.", 0.95);

	public Tree transmissionTree;
	public TraitSet sampleCounts;

	@Override
	public void initAndValidate() {
		transmissionTree = transmissionTreeInput.get();
		sampleCounts = sampleCountsInput.get();
	}

	@Override
	public void run() throws Exception {

		List<Tree> treeList = new ArrayList<>();

		try (PrintStream ps = new PrintStream(fileNameInput.get())) {
			for (int i = 0; i < nSimsInput.get(); i++) {
				Tree tree = new SimulatedGeneTree();
				tree.initByName("transmissionTree", transmissionTreeInput.get(), "sampleCounts",
						sampleCountsInput.get());

				treeList.add(tree);
				ps.println(tree.toString() + ";");
			}
		}

		if (reportFileNameInput.get() != null) {
			try (PrintStream ps = new PrintStream(reportFileNameInput.get())) {
				TreeTraceAnalysis analysis = new TreeTraceAnalysis(treeList, 0.0);
				analysis.computeCredibleSet(credibilityThresholdInput.get());
				analysis.report(ps);
			}
		}

	}


}
