package tnt.util;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;

public class GeneTreeInitializer extends BEASTObject {

	final public Input<Tree> geneTreeInput = new Input<>("geneTree",
			"Gene tree to initialize", Input.Validate.REQUIRED);

	final public Input<TraitSet> sampleCountsInput = new Input<>("sampleCounts",
			"TraitSet defining number of  samples per node in species tree.");

	private Tree geneTree;
	private TraitSet sampleCounts;

	@Override
	public void initAndValidate() {
		geneTree = geneTreeInput.get();
		sampleCounts = sampleCountsInput.get();

	}

	public Tree getGeneTree() {
		return geneTree;
	}

	public TraitSet getSampleCounts() {
		return sampleCounts;
	}

}
