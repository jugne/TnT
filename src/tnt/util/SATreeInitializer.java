package tnt.util;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.RandomTree;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.util.Randomizer;

public class SATreeInitializer extends Tree implements StateNodeInitialiser {

	final public Input<Tree> saTreeInput = new Input<>("saTree",
			"The species tree to initialize.",
			Input.Validate.REQUIRED);

	public Input<List<TaxonSet>> taxonsetsInput = new Input<>("patientTaxonSets",
			"a separate list of taxa for samples collected from the same patient", new ArrayList<>());

	final public Input<RealParameter> birthRate = new Input<>("birthRate",
			"Tree prior birth rate to initialize");

	@Override
	public void initStateNodes() {
		final Tree saTree = saTreeInput.get();
		samePatientSamplingInit(saTree);

	}

	private void samePatientSamplingInit(Tree saTree) {
		final RealParameter birthRateParameter = birthRate.get();
		final Double lambda = (birthRateParameter == null) ? 1.0 : birthRateParameter.getValue();
		final Double initialPopSize = 1.0 / lambda; // scales coalescent tree height inverse to birth rate
		final RealParameter popSize = new RealParameter(initialPopSize.toString());
		final ConstantPopulation pf = new ConstantPopulation();
		pf.setInputValue("popSize", popSize);

		final RandomTree rnd = new RandomTree();
		rnd.setInputValue("taxonset", saTree.getTaxonset());
		if (saTree.hasDateTrait())
			rnd.setInputValue("trait", saTree.getDateTrait());

		rnd.setInputValue("populationModel", pf);
		rnd.setInputValue("populationModel", pf);
		rnd.initAndValidate();

		System.out.println(rnd.getRoot().toNewick());

		/////////////

		int nSets = taxonsetsInput.get().size();
		int nSamples = saTree.getLeafNodeCount();
		List<List<Node>> taxonsets = new ArrayList<List<Node>>(nSets);
		List<Node> singlePatientSamples = new ArrayList<Node>();
		double maxHeight = Double.NEGATIVE_INFINITY;

		// If there are same patient sampling over time supplied
		if (nSets != 0) {
			for (int i = 0; i < nSets; i++)
				taxonsets.add(new ArrayList<Node>());

			for (int idx = 0; idx < nSamples; idx++) {
				Node leaf = rnd.getNode(idx);
				if (leaf.getParent() != null)
					leaf.getParent().removeAllChildren(false);
				int i = 0;
				for (TaxonSet t : taxonsetsInput.get()) {
					List<String> tmp = t.asStringList();
					if (tmp.contains(leaf.getID())) {
						taxonsets.get(i).add(leaf);
						break;
					} else {
						i += 1;
					}
				}
				singlePatientSamples.add(leaf);
				if (leaf.getHeight() > maxHeight)
					maxHeight = leaf.getHeight();
			}

			for (int ii = 0; ii < nSets; ii++) {
				taxonsets.get(ii).sort(Comparator.comparing(Node::getHeight));
			}
		}

		int id = nSamples;
		for (int s = 0; s < nSets; s++) {
			Node[] tmp = new Node[taxonsets.get(s).size()];
			tmp = taxonsets.get(s).toArray(tmp);
			int nNodes = tmp.length;
			for (int n = 0; n < nNodes - 1; n++) {
				Node parent = new Node();
				parent.setID(Integer.toString(id));
				parent.setNr(id);
				id += 1;

				parent.addChild(tmp[n]);
				parent.addChild(tmp[n + 1]);
				singlePatientSamples.remove(tmp[n]);
				singlePatientSamples.remove(tmp[n + 1]);
				parent.setHeight(tmp[n + 1].getHeight());
				tmp[n] = null;
				tmp[n + 1] = parent;
			}
			singlePatientSamples.add(tmp[nNodes - 1]);
			if (tmp[nNodes - 1].getHeight() > maxHeight)
				maxHeight = tmp[nNodes - 1].getHeight();
		}

		while (singlePatientSamples.size() > 1) {
			Node node1 = singlePatientSamples.get(Randomizer.nextInt(singlePatientSamples.size()));
			singlePatientSamples.remove(node1);
			Node node2 = singlePatientSamples.get(Randomizer.nextInt(singlePatientSamples.size()));
			singlePatientSamples.remove(node2);

			double deltaT = Randomizer.nextExponential(
					(singlePatientSamples.size() + 2) * (singlePatientSamples.size() + 1) * 0.5 * 1.0 / initialPopSize);
			double coalTime = maxHeight + deltaT;

			Node parent = new Node();
			parent.setID(Integer.toString(id));
			parent.setNr(id);
			id += 1;
			parent.setHeight(coalTime);
			parent.addChild(node1);
			parent.addChild(node2);

			singlePatientSamples.add(parent);
			maxHeight = coalTime;
		}

		Node root = singlePatientSamples.get(0);

		copyTreeStructure(new Tree(root), saTree);
	}

	// copy the structure of the source tree to the destination tree
	// preserving the leaf node names and numbers
	private void copyTreeStructure(final Tree src, final Tree dst) {
		final Node[] dstNodes = dst.getNodesAsArray();
		final Node[] srcNodes = src.getNodesAsArray();

		final Map<String, Integer> srcTipNumbers = new HashMap<>();

		final int nodeCount = src.getNodeCount();
		final int leafNodeCount = src.getLeafNodeCount();

		// Clear the children of all internal nodes in the destination tree
		for (int nodeNumber = leafNodeCount; nodeNumber < nodeCount; nodeNumber++)
			dstNodes[nodeNumber].removeAllChildren(false);

		// Record the node number of all leaves in the source tree
		for (int nodeNumber = 0; nodeNumber < leafNodeCount; nodeNumber++) {
			final String srcName = srcNodes[nodeNumber].getID();
			srcTipNumbers.put(srcName, nodeNumber);
		}

		// Set the heights of all nodes to match the source height
		for (int nodeNumber = 0; nodeNumber < nodeCount; nodeNumber++) {
			final Node dstNode = dstNodes[nodeNumber];
			Node srcNode;

			// find the corresponding node from the source tree
			if (nodeNumber < leafNodeCount) { // if this is a leaf node
				final String speciesName = dstNode.getID();
				System.out.println(speciesName);
				final int srcTipNumber = srcTipNumbers.get(speciesName);

				srcNode = srcNodes[srcTipNumber];
			} else { // if this is an internal node
				srcNode = srcNodes[nodeNumber];
			}

			// Copy height
			dstNode.setHeight(srcNode.getHeight());

			// Clear and copy metadata
			final Set<String> dstMetaDataNames = dstNode.getMetaDataNames();
			final Set<String> srcMetaDataNames = srcNode.getMetaDataNames();
			final Set<String> srcLengthMetaDataNames = srcNode.getLengthMetaDataNames();

			for (String metaDataName : dstMetaDataNames)
				dstNode.removeMetaData(metaDataName);

			for (String metaDataName : srcMetaDataNames)
				dstNode.setMetaData(metaDataName, srcNode.getMetaData(metaDataName));

			for (String lengthMetaDataName : srcLengthMetaDataNames)
				dstNode.setMetaData(lengthMetaDataName, srcNode.getLengthMetaData(lengthMetaDataName));

			// if this is not the root node, also set the parent and child
			// connections to match the source
			if (nodeNumber != nodeCount - 1) {
				boolean rightChild = false;
				Node siblingNode = null;
				final int parentNumber = srcNode.getParent().getNr();
				final Node parentNode = dstNodes[parentNumber];
				if (srcNode.getParent().getChild(0) == srcNode)
					rightChild = true;
				if (rightChild && parentNode.getChildCount() > 0) {
					siblingNode = parentNode.getChild(0);
					parentNode.removeChild(siblingNode);
				}
				dstNode.setParent(parentNode);
				parentNode.addChild(dstNode);

				if (siblingNode != null) {
					siblingNode.setParent(parentNode);
					parentNode.addChild(siblingNode);
				}

			}
		}
	}

	@Override
	public void getInitialisedStateNodes(final List<StateNode> stateNodes) {
		stateNodes.add(saTreeInput.get());

		final RealParameter brate = birthRate.get();
		if (brate != null) {
			stateNodes.add(brate);
		}
	}

}
