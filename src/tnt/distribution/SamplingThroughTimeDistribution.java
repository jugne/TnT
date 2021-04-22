package tnt.distribution;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import tnt.transmissionTree.TransmissionTree;

public class SamplingThroughTimeDistribution extends Distribution {

	public Input<TransmissionTree> trTreeInput = new Input<>("transmissionTree",
			"transmission tree", Validate.REQUIRED);
	public Input<List<TaxonSet>> taxonsetsInput = new Input<>("taxonsets",
			"a separate list of taxa for samples collected from the same patient");




	private TransmissionTree trTree;
	private List<List<Node>> taxonsets;
	private Integer[] sampleTaxonsetAssignment;
	private int nSets;
	private int nSamples;


	@Override
	public void initAndValidate() {
		trTree = trTreeInput.get();
		nSets = taxonsetsInput.get().size();
		nSamples = trTree.getLeafNodeCount();
		taxonsets = new ArrayList<List<Node>>(nSets);
		sampleTaxonsetAssignment = new Integer[nSamples];

		for (int i = 0; i < nSamples; i++)
			sampleTaxonsetAssignment[i] = -1;

		for (int idx = 0; idx < nSamples; idx++) {
			Node leaf = trTree.getNode(idx);
			int i = 0;
			for (TaxonSet t : taxonsetsInput.get()) {
				List<String> tmp = t.asStringList();
				if (tmp.contains(leaf.getID())) {
					taxonsets.get(i).add(leaf);
					sampleTaxonsetAssignment[leaf.getNr()] = i;
					break;
				} else {
					i += 1;
				}
			}
		}

		for (int ii = 0; ii < nSets; ii++) {
			taxonsets.get(ii).sort(Comparator.comparing(Node::getHeight));
		}

	}



	@Override
	public double calculateLogP() {
		for (int j = 0; j < nSamples; j++) {
			Node startLeaf = trTree.getNode(j);
			int leafNodeAssignment = sampleTaxonsetAssignment[startLeaf.getNr()];

			if (leafNodeAssignment != -1) {
				List<Node> tmp = taxonsets.get(leafNodeAssignment);
				tmp.remove(startLeaf);
				Node child = startLeaf;
				Node parent = startLeaf.getParent();
				while (parent != null) {
					if (parent.isFake()) {
						if (tmp.get(0).getNr() == parent.getChild(1).getNr()) {
							if (!tmp.remove(parent))
								System.out.println("SamplingThroughTime distribution error");
						} else {
							return Double.NEGATIVE_INFINITY;
						}
					} else if (parent.getChild(0) != child) {
						return Double.NEGATIVE_INFINITY;
					}
					child = parent;
					parent = child.getParent();
				}

				if (!tmp.isEmpty())
					return Double.NEGATIVE_INFINITY;
			}
		}

		return 0;
	}


	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub

	}

}
