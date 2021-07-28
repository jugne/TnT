/*
 * Copyright (C) 2019. Tim Vaughan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

package tnt.operators;

import static pitchfork.Pitchforks.getTrueNodes;
import static pitchfork.Pitchforks.isPolytomy;

import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import pitchfork.Pitchforks;
import starbeast2.SpeciesTreeInterface;
import tnt.distribution.GeneTreeIntervals;
import tnt.util.Tools;

@Description("SPR operator for trees with polytomies and multiple mergers.")
// In this version we will assume that all three variants of attachement have equal probability
public class SPROperator extends TreeOperator {

	public Input<Double> rootAttachLambdaInput = new Input<>(
			"rootAttachLambda",
            "Mean of exponential distribution (relative to tree height)" +
                    "used to position attachments above the root.",
			2.0);

	public Input<Double> probBottleneckInput = new Input<>(
			"probBottleneck",
			"Probability of attaching to existing coalescent event or making a multimerger.",
			0.1);

	public Input<SpeciesTreeInterface> transmissionTreeInput = new Input<>("transmissionTree",
			"Fully labeled transmission tree on which to simulate gene trees", Validate.REQUIRED);

	public Input<GeneTreeIntervals> geneTreeIntervalsInput = new Input<>("geneTreeIntervals",
			"intervals for a gene tree", Validate.REQUIRED);

    Tree tree;
	Double rootAttachLambda, probBottleneck;
	Integer[] geneNodeAssignment;
	GeneTreeIntervals intervals;
	SpeciesTreeInterface transmissionTree;

    @Override
    public void initAndValidate() {
		rootAttachLambda = rootAttachLambdaInput.get();
		intervals = geneTreeIntervalsInput.get();
		probBottleneck = probBottleneckInput.get();
		tree = treeInput.get();
		transmissionTree = transmissionTreeInput.get();
    }

    @Override
    public double proposal() {
		double logHR = 0.0;

		boolean bottleneck = false;
		geneNodeAssignment = intervals.getGeneTreeNodeAssignment();

		// Get list of nodes below finite-length edges
		List<Node> trueNodes = getTrueNodes(tree);
		List<Double> trHeights = Tools.getTransmissionHeights(transmissionTree);

		List<Node> rootAndTrueNodesWithParentsAtTransmission = Tools.getGeneRootAndNodesWithParentsAtTransmission(
				trueNodes,
				trHeights);


		// Record number of (true) edges in original tree:
		int nEdges = trueNodes.size() - rootAndTrueNodesWithParentsAtTransmission.size();

		// Select non-root subtree at random

		Node srcNode, srcNodeParent;
		do {
			srcNode = trueNodes.get(Randomizer.nextInt(trueNodes.size()));
			srcNodeParent = srcNode.getParent();
		} while (srcNodeParent == null || rootAndTrueNodesWithParentsAtTransmission.contains(srcNode)); // cannot detach
																										// nodes
																									// at
																							// transmission event,
																							// because this operator
																							// cannot put them there.
																							// There is
		// TransmissionAttach
																							// operator for that.

		int trNodeNr = geneNodeAssignment[srcNode.getNr()];
		List<Integer> possibleAssignments = new ArrayList<>();
		possibleAssignments.add(trNodeNr);
		possibleAssignments.addAll(Tools.getAllParentNrs(transmissionTree.getNode(trNodeNr)));
		possibleAssignments.addAll(Tools.getChildNrs(transmissionTree.getNode(trNodeNr)));

		Node logicalParent = Pitchforks.getLogicalParent(srcNode);
		Node srcNodeSister = getOtherChild(srcNodeParent, srcNode);

		// Record whether the the original attachment was a polytomy or merger
		boolean origAttachWasPolytomy = isPolytomy(srcNodeParent);
		boolean origAttachWasMerger = Tools.isMultiMerger(trueNodes, logicalParent);
		boolean origAttachWasBottleneck = origAttachWasPolytomy || origAttachWasMerger;

		boolean parentWasRoot = srcNodeParent.isRoot();
		double oldParentHeight = srcNodeParent.getHeight();

		// Disconnect subtree

		srcNodeParent.removeChild(srcNodeSister);

		Node srcNodeGrandparent = null;
		if (parentWasRoot) {
			srcNodeSister.setParent(null);
		} else {
			srcNodeGrandparent = srcNodeParent.getParent();
			srcNodeGrandparent.removeChild(srcNodeParent);
			srcNodeGrandparent.addChild(srcNodeSister);
		}

		srcNodeParent.setParent(null);

		// Select new attachment node

		Node remainingSubtreeRoot;
		if (parentWasRoot)
			remainingSubtreeRoot = srcNodeSister;
		else
			remainingSubtreeRoot = tree.getRoot();

		List<Node> subtreeNodes = getFitNodesInSubtree(remainingSubtreeRoot, srcNode.getHeight(), possibleAssignments);
		if (subtreeNodes.size() == 0)
			return Double.NEGATIVE_INFINITY; // no nodes to choose from
//		List<Node> innerNodes = Pitchforks.getTrueInternalNodes(tree);

		List<Node> subtreeNodesAtTransmission = Tools.getGeneNodesAtTransmissionWithPrecision(subtreeNodes,
				trHeights);
		
//		if (Tools.listEqualsIgnoreOrder(subtreeNodes, subtreeNodesAtTransmission)) {
//			return Double.NEGATIVE_INFINITY;
//		}

		List<Node> subtreeNodesBottlenecHeight = new ArrayList<Node>();
		for (Node n : subtreeNodes) {
			if (!subtreeNodesAtTransmission.contains(n) && Tools.greaterDouble(n.getHeight(), srcNode.getHeight())
					&& !n.isLeaf()) {
				subtreeNodesBottlenecHeight.add(n);
			}
		}

		Node heightNode;
		heightNode = subtreeNodes.get(Randomizer.nextInt(subtreeNodes.size()));



		Node newAttachNode;
		Double newHeight = null;
		Node heightNodeBot = null;

		List<Node> fitNodes = new ArrayList<Node>();
		List<Node> atSameHeight = new ArrayList<Node>();

		if (subtreeNodesBottlenecHeight.size() > 0) {
			heightNodeBot = subtreeNodesBottlenecHeight
					.get(Randomizer.nextInt(subtreeNodesBottlenecHeight.size()));
			newHeight = heightNodeBot.getHeight();

			for (Node n : subtreeNodes) {
				if (!subtreeNodesAtTransmission.contains(n)) { // do not allow exact height
			// to be transmission height
					if (Tools.equalWithPrecisionDouble(n.getHeight(), newHeight)
							&& subtreeNodesBottlenecHeight.contains(n))
					atSameHeight.add(n);
					if (n.getHeight() <= newHeight
							&& (n.isRoot() || Pitchforks.getLogicalParent(n).getHeight() > newHeight)) {
					fitNodes.add(n);
				}
		}
			}
		} else {
			bottleneck = false;
		}

		if (fitNodes.size() == 0 || heightNodeBot.isLeaf()
				|| Tools.greaterOrEqualDouble(srcNode.getHeight(), newHeight))
			bottleneck = false;
		else {
			if (Randomizer.nextDouble() < probBottleneck) {
				bottleneck = true;
				logHR -= Math.log(probBottleneck);
			} else {
				logHR -= Math.log(1 - probBottleneck);

			}
		}

		if (bottleneck) {
			// select a new attachment node
			newAttachNode = fitNodes.get(Randomizer.nextInt(fitNodes.size()));
			// account for its probability
			logHR -= Math.log(1.0 / fitNodes.size());
			// account for probability to pick this specific height:
			// nodes at the same height / all possible nodes
			logHR -= Math.log(atSameHeight.size() / (double) subtreeNodesBottlenecHeight.size()); // (subtreeNodes.size()
																											// -
																											// subtreeNodesAtTransmission.size()));
		} else {
			// if attaching at uniform height, we do not care about transmission heights,
			// since the probability of picking them is vanishing
			logHR -= Math.log(1.0 / (subtreeNodes.size()));
			newAttachNode = heightNode;

			// pick new height and account for its probability
			if (newAttachNode.isRoot()) {
				double minAttachheight = Tools.findMinHeight(newAttachNode, possibleAssignments, geneNodeAssignment,
						transmissionTree);
				double offset = Math.max(minAttachheight, srcNode.getHeight());
//				double offset = Math.max(srcNode.getHeight(), newAttachNode.getHeight());
				double expRate = 1.0 / (rootAttachLambda * offset);
				newHeight = offset + Randomizer.nextExponential(expRate);

				logHR -= -expRate * (newHeight - offset)
						+ Math.log(expRate);
			} else {
				double sibHeight = Tools.findMinHeight(newAttachNode, possibleAssignments, geneNodeAssignment,
						transmissionTree);
				;
				double L = newAttachNode.getParent().getHeight() -
						Math.max(srcNode.getHeight(), sibHeight);
				newHeight = Randomizer.nextDouble() * L +
						Math.max(srcNode.getHeight(), sibHeight);

				logHR -= Math.log(1.0 / L);
			}
		}

		// Incorporate probability of existing attachment point into HR
		List<Node> origFitNodes = new ArrayList<Node>();
		List<Node> origAtSameHeight = new ArrayList<Node>();
		for (Node n : subtreeNodes) {
			if (!subtreeNodesAtTransmission.contains(n)) {
				if (Tools.equalWithPrecisionDouble(n.getHeight(), oldParentHeight)
						&& subtreeNodesBottlenecHeight.contains(n))
					origAtSameHeight.add(n);
				if (Tools.greaterOrEqualDouble(oldParentHeight, n.getHeight())
						&& (n.isRoot()
								|| Tools.greaterDouble(Pitchforks.getLogicalParent(n).getHeight(), oldParentHeight))) {
					origFitNodes.add(n);
				}
			}
		}

		if (origAttachWasBottleneck) {
			logHR += Math.log(probBottleneck);
			logHR += Math.log(1.0 / origFitNodes.size());
			logHR += Math
					.log(origAtSameHeight.size() / (double) subtreeNodesBottlenecHeight.size()); // (subtreeNodes.size()
																											// -
																											// subtreeNodesAtTransmission.size()));
		} else {
			if ((origFitNodes.size() == 1 && !srcNodeSister.isLeaf())
					|| (origFitNodes.size() > 1) && oldParentHeight > srcNode.getHeight()) {
				logHR += Math.log(1 - probBottleneck);
			} else {

//				if (subtreeNodes.contains(srcNodeSister)) {
//					System.out.println();
//				}
			}

			logHR += Math.log(1.0 / subtreeNodes.size());

			if (parentWasRoot) {
				double minAttachheight = Tools.findMinHeight(srcNodeSister, possibleAssignments, geneNodeAssignment,
						transmissionTree);
				double offset = Math.max(minAttachheight, srcNode.getHeight());
//				double offset = Math.max(srcNodeSister.getHeight(), srcNode.getHeight());
				double expRate = 1.0 / (rootAttachLambda * offset);
				logHR += -expRate * (oldParentHeight - offset) + Math.log(expRate);
			} else {
				double sibHeight = Tools.findMinHeight(srcNodeSister, possibleAssignments, geneNodeAssignment,
						transmissionTree);
				double L = srcNodeGrandparent.getHeight()
						- Math.max(sibHeight, srcNode.getHeight());
				logHR += Math.log(1.0 / L);
			}
		}

		// Reconnect subtree

		srcNodeParent.setHeight(newHeight);

		if (newAttachNode.isRoot()) {
			srcNodeParent.addChild(newAttachNode);
		} else {
			Node oldParent = newAttachNode.getParent();
			oldParent.removeChild(newAttachNode);
			oldParent.addChild(srcNodeParent);
			srcNodeParent.addChild(newAttachNode);
		}

		// Ensure correct root if set if this has been modified:
		if (srcNodeSister.isRoot())
			tree.setRoot(srcNodeSister);
		else if (srcNodeParent.isRoot())
			tree.setRoot(srcNodeParent);

//		// Reassign nodes and check compatibility after the move
//		geneNodeAssignment.clear();
//		geneNodeAssignment.putAll(intervals.geneTreeTipAssignment);
//
//		if (!Tools.fillAssignmentAndCheck(intervals.transmissionTreeInput.get(), tree.getRoot(),
//				geneNodeAssignment, null))
//			return Double.NEGATIVE_INFINITY; // not compatible


		List<Node> trueNodesAfter = Pitchforks.getTrueNodes(tree);
		List<Node> rootAndTrueNodesWithParentsAtTransmissionAfter = Tools.getGeneRootAndNodesWithParentsAtTransmission(
				trueNodesAfter,
				trHeights);

		int nEdgesAfter = trueNodesAfter.size() - rootAndTrueNodesWithParentsAtTransmissionAfter.size();


		logHR -= Math.log(1.0 / nEdges);
		logHR += Math.log(1.0 / nEdgesAfter);



		return logHR;
    }

	/**
	 * Gets nodes with parent edges above or including minAge time of the tree.
	 * 
	 * @param subtreeRoot root node at which start looking
	 * @param minAge      minimum possible attachment age
	 * @return list of nodes, which have parent edges above or including minAge.
	 */
	private List<Node> getFitNodesInSubtree(Node subtreeRoot, double minAge, List<Integer> possibleAssignments) {
		List<Node> nodeList = new ArrayList<>();

		if (subtreeRoot.isRoot() || (subtreeRoot.getParent().getHeight() > subtreeRoot.getHeight())
				&& possibleAssignments.contains(geneNodeAssignment[subtreeRoot.getParent().getNr()]))
			nodeList.add(subtreeRoot);

		if (subtreeRoot.getHeight() > minAge) {
			for (Node child : subtreeRoot.getChildren())
				nodeList.addAll(getFitNodesInSubtree(child, minAge, possibleAssignments));
		}

		return nodeList;
	}

//	private List<Node> getGeneNodesAtTransmission(List<Node> nodeList, List<Double> transmissionHeights) {
//		List<Node> tmp = new ArrayList<>();
//		for (Node n : nodeList) {
//			if (transmissionHeights.contains(n.getHeight()))
//				tmp.add(n);
//		}
//
//		return tmp;
//	}

	private List<Double> getHeightsForMultiMerger(List<Node> nodesInSubtree, Node newAttach, double minAge) {
		List<Double> heights = new ArrayList<Double>();

		if (newAttach.isRoot())
			return heights;
		for (Node n : nodesInSubtree) {
			if (!n.isLeaf() && n.getNr() != newAttach.getNr()
					&& n.getHeight() < Pitchforks.getLogicalParent(newAttach).getHeight()
					&& n.getHeight() > newAttach.getHeight() && n.getHeight() > minAge) {
				heights.add(n.getHeight());
			}
		}
		return heights;
	}





}
