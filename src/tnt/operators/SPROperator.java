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

    Tree tree;
	Double rootAttachLambda, probBottleneck;
	SpeciesTreeInterface transmissionTree;

    @Override
    public void initAndValidate() {
		rootAttachLambda = rootAttachLambdaInput.get();
		probBottleneck = probBottleneckInput.get();
		tree = treeInput.get();
		transmissionTree = transmissionTreeInput.get();
    }

    @Override
    public double proposal() {
		double logHR = 0.0;

		boolean bottleneck = false;

		// Get list of nodes below finite-length edges
		List<Node> trueNodes = getTrueNodes(tree);
		List<Double> trHeights = Tools.getTransmissionHeights(transmissionTree);

		List<Node> trueNodesWithParentsAtTransmission = Tools.getGeneNodesWithParentsAtTransmission(trueNodes,
				trHeights);

		// Record number of (true) edges in original tree:
		int nEdges = trueNodes.size() - 1 - trueNodesWithParentsAtTransmission.size();

		// Select non-root subtree at random

		Node srcNode, srcNodeParent;
		do {
			srcNode = trueNodes.get(Randomizer.nextInt(trueNodes.size()));
			srcNodeParent = srcNode.getParent();
		} while (srcNodeParent == null || trueNodesWithParentsAtTransmission.contains(srcNode)); // cannot detach nodes
																									// at
																							// transmission event,
																							// because this operator
																							// cannot put them there.
																							// There is
																							// TRansmissionAttach
																							// operator for that.

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
		if (srcNodeParent.isRoot()) {
			srcNodeSister.setParent(null);
		} else {
			srcNodeGrandparent = srcNodeParent.getParent();
			srcNodeGrandparent.removeChild(srcNodeParent);
			srcNodeGrandparent.addChild(srcNodeSister);
		}

		srcNodeParent.setParent(null);

		// Select new attachment node

		Node remainingSubtreeRoot;
		if (srcNodeSister.isRoot())
			remainingSubtreeRoot = srcNodeSister;
		else
			remainingSubtreeRoot = tree.getRoot();

		List<Node> subtreeNodes = getNodesInSubtree(remainingSubtreeRoot, srcNode.getHeight());
		List<Node> innerNodes = Pitchforks.getTrueInternalNodes(tree);

		List<Node> subtreeNodesAtTransmission = Tools.getGeneNodesAtTransmission(subtreeNodes,
				trHeights);
		
		if (Tools.listEqualsIgnoreOrder(subtreeNodes, subtreeNodesAtTransmission)) {
			return Double.NEGATIVE_INFINITY;
		}

		Node heightNode;
		heightNode = subtreeNodes.get(Randomizer.nextInt(subtreeNodes.size()));


//		while (heightNode.getHeight() == srcNode.getHeight()) {
//			heightNode = subtreeNodes.get(Randomizer.nextInt(subtreeNodes.size()));
//		}
		List<Node> eligibleSubtreeNodes = new ArrayList<>();
		for (Node node : subtreeNodes) {
			if (!node.isLeaf() || node.getHeight() > srcNode.getHeight())
				eligibleSubtreeNodes.add(node);
		}
//
		Node newAttachNode;
		Double newHeight;

			List<Node> fitNodes = new ArrayList<Node>();
		List<Node> atSameHeight = new ArrayList<Node>();
			newHeight = heightNode.getHeight();
		for (Node n : subtreeNodes) {
			if (!subtreeNodesAtTransmission.contains(n)) { // do not allow exact height to be transmission height
				if (n.getHeight() == newHeight)
					atSameHeight.add(n);
				if (n.getHeight() <= newHeight
						&& (n.isRoot() || Pitchforks.getLogicalParent(n).getHeight() > newHeight)) {
					fitNodes.add(n);
				}
			}
		}
		if (fitNodes.size() == 0 || heightNode.isLeaf() || newHeight <= srcNode.getHeight())
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
			logHR -= Math.log(atSameHeight.size() / (double) (subtreeNodes.size() - subtreeNodesAtTransmission.size()));
		} else {
			// if attaching at uniform height, we do not care about transmission heights,
			// since the probability of picking them is vanishing
			logHR -= Math.log(1.0 / (subtreeNodes.size()));
			newAttachNode = heightNode;

			// pick new height and account for its probability
			if (newAttachNode.isRoot()) {
				double offset = Math.max(srcNode.getHeight(), newAttachNode.getHeight());
				double expRate = 1.0 / (rootAttachLambda * offset);
				newHeight = offset + Randomizer.nextExponential(expRate);

				logHR -= -expRate * (newHeight - offset)
						+ Math.log(expRate);
			} else {
				double L = newAttachNode.getParent().getHeight() -
						Math.max(srcNode.getHeight(), newAttachNode.getHeight());
				newHeight = Randomizer.nextDouble() * L +
						Math.max(srcNode.getHeight(), newAttachNode.getHeight());

				logHR -= Math.log(1.0 / L);
			}
		}

		// Incorporate probability of existing attachment point into HR
		List<Node> origFitNodes = new ArrayList<Node>();
		List<Node> origAtSameHeight = new ArrayList<Node>();
		for (Node n : subtreeNodes) {
			if (!subtreeNodesAtTransmission.contains(n)) {
				if (n.getHeight() == oldParentHeight)
					origAtSameHeight.add(n);
				if (n.getHeight() <= oldParentHeight
					&& (n.isRoot() || Pitchforks.getLogicalParent(n).getHeight() > oldParentHeight)) {
					origFitNodes.add(n);
				}
			}
		}

		if (origAttachWasBottleneck) {
			logHR += Math.log(probBottleneck);
			logHR += Math.log(1.0 / origFitNodes.size());
			logHR += Math
					.log(origAtSameHeight.size() / (double) (subtreeNodes.size() - subtreeNodesAtTransmission.size()));
		} else {
			if ((origFitNodes.size() == 1 && !srcNodeSister.isLeaf())
					|| (origFitNodes.size() > 1) && oldParentHeight > srcNode.getHeight()) {
				logHR += Math.log(1 - probBottleneck);

			}
			logHR += Math.log(1.0 / subtreeNodes.size());

			if (parentWasRoot) {
				double offset = Math.max(srcNodeSister.getHeight(), srcNode.getHeight());
				double expRate = 1.0 / (rootAttachLambda * offset);
				logHR += -expRate * (oldParentHeight - offset) + Math.log(expRate);
			} else {
				double L = srcNodeGrandparent.getHeight()
						- Math.max(srcNodeSister.getHeight(), srcNode.getHeight());
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


		List<Node> trueNodesAfter = getTrueNodes(tree);
		List<Node> trueNodesWithParentsAtTransmissionAfter = Tools.getGeneNodesWithParentsAtTransmission(trueNodesAfter,
				trHeights);

		int nEdgesAfter = trueNodesAfter.size() - 1 - trueNodesWithParentsAtTransmissionAfter.size();

		logHR -= Math.log(1.0 / nEdges);
		logHR += Math.log(1.0 / nEdgesAfter);


		return logHR;
    }

	private List<Node> getNodesInSubtree(Node subtreeRoot, double minAge) {
		List<Node> nodeList = new ArrayList<>();

		if (subtreeRoot.isRoot() || subtreeRoot.getParent().getHeight() > subtreeRoot.getHeight())
			nodeList.add(subtreeRoot);

		if (subtreeRoot.getHeight() > minAge) {
			for (Node child : subtreeRoot.getChildren())
				nodeList.addAll(getNodesInSubtree(child, minAge));
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
