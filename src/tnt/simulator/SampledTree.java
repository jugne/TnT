package tnt.simulator;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import bdmmprime.parameterization.Parameterization;
import beast.core.Function;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import tnt.util.Tools;

//TODO rho sampling, r>0
public class SampledTree extends Tree {
	public Input<Tree> treeInput = new Input<>("tree", "Full Tree",
			Input.Validate.REQUIRED);

	public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
			"BDMM parameterization",
			Input.Validate.REQUIRED);

	public Input<Function> finalSampleOffsetInput = new Input<>("finalSampleOffset",
			"If provided, the difference in time between the final sample and the end of the BD process.",
			new RealParameter("0.0"));

	public Input<Function> hiddenEventsCounterInput = new Input<>("nHiddenEvents",
			"Number of hidden events. Set by sampler.",
			new RealParameter("0.0"));

	private RealParameter finalSampleOffsetFullTree;
	public Function finalSampleOffsetSampledTree;
	private RealParameter hiddenEventsCounter;

	private Tree fullTree;
	private double psi_i;
	private double r_i;
	private double rho_i;
	private double processLength;
	private Parameterization parameterization;

	private double youngestLeafAge;
	private HashSet<Node> sampledNodes;
	private int sampleCount;
	private Double nHiddenEvents;

	HashMap<Integer, Integer> hostSamples;

	private int nr;

	public void initAndValidate() {
		sampleCount = 0;
		nHiddenEvents = 0.;
		sampledNodes = new HashSet<Node>();
		hostSamples = new HashMap<Integer, Integer>();
		youngestLeafAge = Double.MAX_VALUE;
		parameterization = parameterizationInput.get();
		processLength = parameterization.getTotalProcessLength();
		finalSampleOffsetFullTree = (RealParameter) finalSampleOffsetInput.get();
		hiddenEventsCounter = (RealParameter) hiddenEventsCounterInput.get();

//	}
//
//	public void sample() {
//		origin = parameterization.originInput.get().getArrayValue(0);
		fullTree = treeInput.get();
//		finalSampleOffsetFullTree = finalSampleOffsetInput.get();

		hostSamples.put(Integer.valueOf(fullTree.getRoot().getID().split("_")[0]), 0);
		nr = fullTree.getRoot().getNodeCount();
		Node sampledRoot = sampleTree(fullTree.getRoot());
		// traverse the tree from sampled nodes towards the root and mark
		// the nodes that belong to the sampled tree
		HashSet<Node> parents;
		HashSet<Node> children = sampledNodes;
		while (children.size() > 1) {
			parents = collectParents(children);
			children = parents;
		}
		removeSingleChildNodes(sampledRoot);
		if (sampledRoot.getChildCount() == 1 && sampledRoot.getNr() != 0) {
			Node newRoot = sampledRoot.getLeft();
			sampledRoot = newRoot;
		}
		Tools.numberInternalNodesOnSubtree(sampledRoot, sampledNodes.size());
//		Integer i = 0;
//		for (Node l : sampledRoot.getAllLeafNodes()) {
//			l.setID(i.toString());
//			i++;
//		}

		assignFromWithoutID(new Tree(sampledRoot));

		List<Node> sampledTreeLeafNodes = sampledRoot.getAllLeafNodes();
		List<String> sampledTreeLeafNodeIds = new ArrayList<String>();

		for (Node k : sampledTreeLeafNodes) {
			sampledTreeLeafNodeIds.add(k.getID().split("_")[0]);

		}
		updateHiddenEventsCounter(fullTree.getRoot(), sampledTreeLeafNodeIds);

		finalSampleOffsetFullTree.setValue(youngestLeafAge);
		hiddenEventsCounter.setValue(nHiddenEvents);

	}

	private Node sampleTree(Node subroot) {
		Node root = null;
		Node currentNode = new Node();
		currentNode.setHeight(subroot.getHeight());
		currentNode.setID(subroot.getID());
//		if (subroot.isLeaf()) {
		currentNode.setNr(-1);
//		}

		double t_start_branch = processLength;

		if (!subroot.isRoot())
			t_start_branch = subroot.getParent().getHeight();
		double t_end_branch = subroot.getHeight();

		int i = parameterization.getNodeIntervalIndex(subroot, finalSampleOffsetFullTree.getArrayValue());
		double t_end_int = 0;
		if (parameterization.getTotalIntervalCount() != 1)
			t_end_int = parameterization.getAge(parameterization.getIntervalEndTimes()[i - 1],
					finalSampleOffsetFullTree.getArrayValue());

		updateParameters(i);

		boolean first = true;
		double t_start = t_start_branch;
		while (t_end_branch < t_end_int) {
			// record event times
			double deltaT = Randomizer.nextExponential(psi_i);
			boolean stop = false;
			while (t_end_int < t_start - deltaT && t_end_branch < t_start - deltaT) {
				// need to fix
				Node tmp = new Node();
				if (!currentNode.isRoot()) {
					Node previousParent = currentNode.getParent();
					Tools.replaceNodeKeepDirection(previousParent, currentNode, tmp);
//					tmp.setID(previousParent.getID());
//					previousParent.removeChild(currentNode);
//					previousParent.addChild(tmp);
				}
				tmp.setHeight(t_start - deltaT);
				if (first) {
					root = tmp;
					first = false;
				}
				tmp.setNr(sampledNodes.size());
				sampleCount++;
				sampledNodes.add(tmp);
				if (Randomizer.nextDouble() < r_i) {
					stop = true;
					break;
				}

				tmp.addChild(currentNode);

				t_start -= deltaT;

				deltaT = Randomizer.nextExponential(psi_i);

			}
			if (stop)
				break;

			i = i - 1;
			t_start = t_end_int;
			t_end_int = parameterization.getIntervalEndTimes()[i - 1];
			updateParameters(i);
		}

		double deltaT = Randomizer.nextExponential(psi_i);
		if (first && t_end_branch > t_start - deltaT) {
			root = currentNode;
		}
		while (t_end_branch < t_start - deltaT) {
			Node tmp = new Node();
			String tmpId = currentNode.getID().split("_")[0];
			if (!currentNode.isRoot()) {
				Node previousParent = currentNode.getParent();
				Tools.replaceNodeKeepDirection(previousParent, currentNode, tmp);
				tmpId = previousParent.getID().split("_")[0];
//				previousParent.removeChild(currentNode);
//				previousParent.addChild(tmp);
			}
			tmp.setHeight(t_start - deltaT);
			if (first) {
				root = tmp;
				first = false;
			}


			if (Randomizer.nextDouble() > r_i) {
				Node right = new Node();
				right.setNr(sampleCount);
				sampleCount++;
				sampledNodes.add(right);
				right.setHeight(t_start - deltaT);
				tmp.setLeft(currentNode);
				tmp.setRight(right);
				currentNode.setParent(tmp);
				right.setParent(tmp);

				if (right.getHeight() < youngestLeafAge) {
					youngestLeafAge = right.getHeight();
				}

				Integer tmpNr = hostSamples.get(Integer.valueOf(tmpId));
				if (tmpNr != null)
					tmpNr += 1;
				else
					tmpNr = 1;
//				int tmpNr = hostSamples.get(Integer.valueOf(tmpId)) + 1;
				right.setID(tmpId + "_" + tmpNr + "_");
				tmp.setID(tmpId + "_");
				tmp.setNr(nr);
				nr++;
				hostSamples.put(Integer.valueOf(tmpId), tmpNr);

			} else {
				tmp.setNr(sampledNodes.size());
				sampleCount++;
				sampledNodes.add(tmp);
				if (tmp.getHeight() < youngestLeafAge) {
					youngestLeafAge = tmp.getHeight();
				}
				int tmpNr = hostSamples.get(Integer.valueOf(tmpId)) + 1;
				tmp.setID(tmpId + "_" + tmpNr + "_");
				hostSamples.put(Integer.valueOf(tmpId), tmpNr);

				return root;
//				tmp.addChild(currentNode);
			}


			t_start -= deltaT;

			deltaT = Randomizer.nextExponential(psi_i);

		}

		if (!subroot.isLeaf()) {
			currentNode.addChild(sampleTree(subroot.getChild(0)));
			currentNode.addChild(sampleTree(subroot.getChild(1)));
		}
//		else if (subroot.getHeight() < youngestLeafAge) {
//			youngestLeafAge = subroot.getHeight();
//		}

//		root = currentNode;
		return root;
	}

	/**
	 * collect parents of nodes in children set. During the collection all visited
	 * nodes get non-negative numbers in order to extract the sampled tree later.
	 * 
	 * @param children set of nodes
	 * @return parents of nodes in children set
	 */
	public HashSet<Node> collectParents(HashSet<Node> children) {
		HashSet<Node> parents = new HashSet<Node>();
		for (Node node : children) {
			if (node.getParent() != null) {
				if (node.getParent().getNr() == -1) {
					node.getParent().setNr(sampleCount);
					sampleCount++;
				}
				parents.add(node.getParent());
			} else {
				parents.add(node);
			}
		}
		return parents;
	}

	/**
	 * Extract the sampled tree by discarding all the nodes that have -1 number
	 * simultaneously suppress single child nodes (nodes with only one child
	 * numbered by non-negative number)
	 * 
	 * @param node
	 */
	public void removeSingleChildNodes(Node node) {
		if (node != null && !node.isLeaf()) {
			Node left = node.getLeft();
			Node right = node.getRight();

			removeSingleChildNodes(right);
			removeSingleChildNodes(left);

			List<Node> toRemove = new ArrayList<>();

			for (Node child : node.getChildren()) {
				if (child.getNr() == -1) {
					toRemove.add(child);
//					node.removeChild(child);
				}
			}
			for (Node n : toRemove) {
				node.removeChild(n);
			}

			if (node.getChildCount() == 1 && node.getParent() != null) {
				Node parent = node.getParent();
				Node newChild = node.getLeft();
				Boolean nodeIsLeft = node.getNr() == parent.getLeft().getNr() ? true : false;
				Node otherSibling = nodeIsLeft ? parent.getRight() : parent.getLeft();
				parent.removeChild(node);
				parent.removeChild(otherSibling);
				if (nodeIsLeft) {
					parent.addChild(newChild);
					parent.addChild(otherSibling);
				} else {
					parent.addChild(otherSibling);
					parent.addChild(newChild);
				}
				newChild.setParent(parent);
				otherSibling.setParent(parent);
				// node.setParent(null);
			}
		}
	}

	private void updateParameters(int intervalNr) {
		psi_i = parameterization.getSamplingRates()[intervalNr][0];
		rho_i = parameterization.getRhoValues()[intervalNr][0];
		r_i = parameterization.getRemovalProbs()[intervalNr][0];
	}
	
	private void updateHiddenEventsCounter(Node fullTreeSubroot, List<String> sampledTreeLeafNodeIds) {
		if (!fullTreeSubroot.isLeaf()) {
			boolean hidden = false;
			if (!fullTreeSubroot.isFake()) {
				List<Node> fullTreeLeavesRight = new ArrayList<>();
				fullTreeSubroot.getRight().getAllLeafNodes(fullTreeLeavesRight);
				List<Node> fullTreeLeavesLeft = new ArrayList<>();
				fullTreeSubroot.getLeft().getAllLeafNodes(fullTreeLeavesLeft);

				for (Node n : fullTreeLeavesRight) {
					if (sampledTreeLeafNodeIds.contains(n.getID().split("_")[0])) {
						hidden = true;
						break;
					}
				}

				if (hidden) {
					for (Node n : fullTreeLeavesLeft) {
						if (sampledTreeLeafNodeIds.contains(n.getID().split("_")[0])) {
							hidden = false;
							break;
						}
//						hidden = false;
					}
				}
				if (hidden) {
					nHiddenEvents += 1;
//					writerHiddenEventsTime.print("\t" + (rhoSamplingTime
//							+ fullTreeSubroot.getHeight()));
				}

			}
			updateHiddenEventsCounter(fullTreeSubroot.getRight(), sampledTreeLeafNodeIds);
			updateHiddenEventsCounter(fullTreeSubroot.getLeft(), sampledTreeLeafNodeIds);

		}
	}

	@Override
	public void init(PrintStream out) {
//        untypedTree.init(out);
	}

    public void log(long sample, PrintStream out) {
        Tree tree = (Tree) getCurrent();
//        out.print("tree STATE_" + sample + " = ");
        final int[] dummy = new int[1];
		final String newick = tree.getRoot().toNewick();
        out.print(newick);
        out.print(";");
    }

}
