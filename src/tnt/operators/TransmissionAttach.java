/*
 * Copyright (C) 2020. Ugne Stolz
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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import pitchfork.Pitchforks;
import tnt.distribution.GeneTreeIntervals;
import tnt.util.Tools;

public class TransmissionAttach extends TreeOperator {

	public Input<GeneTreeIntervals> geneTreeIntervalsInput = new Input<>("geneTreeIntervals",
			"intervals for a gene tree", Validate.REQUIRED);

	public Input<Double> rootAttachLambdaInput = new Input<>(
			"rootAttachLambda",
			"Mean of exponential distribution (relative to tree height)" +
					"used to position attachments above the root.",
			2.0);

	Tree geneTree;
	GeneTreeIntervals intervals;
	Integer[] geneNodeAssignment;
	HashMap<Integer, List<Integer>> inverseGeneNodeAssignment;
	double[] trNodeOccupancy;
	Double rootAttachLambda;
	double[] trHeights;

	int nGeneNodes;
	int nTrNodes;

    @Override
    public void initAndValidate() {
		geneTree = treeInput.get();
		intervals = geneTreeIntervalsInput.get();
		rootAttachLambda = rootAttachLambdaInput.get();

		nGeneNodes = geneTree.getNodeCount();
		nTrNodes = intervals.transmissionTreeInput.get().getNodeCount();

		trNodeOccupancy = new double[nGeneNodes * nTrNodes];
		geneNodeAssignment = new Integer[geneTree.getNodeCount()];

    }

    @Override
    public double proposal() {

        double logHR = 0.0;
		System.arraycopy(intervals.getGeneTreeNodeAssignment(), 0, geneNodeAssignment, 0,
				intervals.getGeneTreeNodeAssignment().length);

//		inverseGeneNodeAssignment = new HashMap<Integer, List<Integer>>(intervals.getInverseGeneTreeNodeAssignment());
		// fill the recipient node list on
		List<Node> recipients = new ArrayList<>();
		Node trTreeRoot = intervals.transmissionTreeInput.get().getRoot();
		getRecipients(trTreeRoot, recipients);

		// choose a transmission tree recipient node
		// transmission event happens at its parent node
		int nRecipients = recipients.size();
		if (nRecipients == 0)
			return Double.NEGATIVE_INFINITY;
		Node chosenRecipient = recipients.get(Randomizer.nextInt(nRecipients));
		// store numbers of all nodes that are descending from the chosen recipient node
		List<Integer> recipientGroupNr = new ArrayList<>();
		recipientGroupNr.add(chosenRecipient.getNr());
		recipientGroupNr.addAll(getChildNrs(chosenRecipient));

		List<Integer> possibleAssignments = new ArrayList<>();
		possibleAssignments.addAll(getAllParentNrs(chosenRecipient));
		possibleAssignments.addAll(recipientGroupNr);

		List<Node> trueNodes = Pitchforks.getTrueNodes(geneTree);
		trHeights = intervals.getTrHeights();
//		List<Node> trueNodesWithParentsNotAtTransmission = Tools.getGeneNodesWithParentNotAtTransmission(trueNodes,
//				trHeights);
		List<Node> fitToMove = new ArrayList<>();
		List<Node> fitToMoveAfter = new ArrayList<>();
		List<Node> fitToAttach = new ArrayList<>();
		List<Node> fitToAttachAfter = new ArrayList<>();

		if (Randomizer.nextBoolean()) {
			// attach a node at transmission
//			System.out.println("attach at: " + chosenRecipient.getParent().getHeight());
			

			// get all nodes, that are not at transmission height and can be attached there
			fitToMove = getNodesToMoveToTransmission(trueNodes,
					chosenRecipient,
					recipientGroupNr);
			if (fitToMove.size() == 0)
				return Double.NEGATIVE_INFINITY;
			
			// Choose node which will be moved, store its parent and sibling nodes;
			// parent node is the one we will detach and re-attach.
			Node srcNode = fitToMove.get(Randomizer.nextInt(fitToMove.size()));
			Node srcNodeParent = srcNode.getParent();
			Node srcNodeSibling = getOtherChild(srcNode.getParent(), srcNode);
			Node srcNodeGrandparent = null;

			// Account for current position probability

			if (srcNodeParent.isRoot()) {
//				double offset;
				double minAttachheight = findMinHeight(srcNodeSibling, possibleAssignments);
				double offset = Math.max(minAttachheight, srcNode.getHeight());
//				if (possibleLocations.contains(geneNodeAssignment.get(srcNodeSibling.getNr())))
//					offset = Math.max(srcNodeSibling.getHeight(), srcNode.getHeight());
//				else {
//					int parentAssignmentLabel = geneNodeAssignment.get(srcNodeSibling.getParent().getNr());
//					Node parentAssigntmentNode = intervals.transmissionTreeInput.get().getNode(parentAssignmentLabel);
//					offset = Math.max(parentAssigntmentNode.getHeight(), srcNode.getHeight());
//				}
				double expRate = 1.0 / (rootAttachLambda * offset);
				logHR += -expRate * (srcNodeParent.getHeight() - offset) + Math.log(expRate);
			} else {
				double sibHeight = findMinHeight(srcNodeSibling, possibleAssignments);
//				if (possibleLocations.contains(geneNodeAssignment.get(srcNodeSibling.getNr())))
//					sibHeight = srcNodeSibling.getHeight();
//				else {
//					int parentAssignmentLabel = geneNodeAssignment.get(srcNodeSibling.getParent().getNr());
//					Node parentAssigntmentNode = intervals.transmissionTreeInput.get().getNode(parentAssignmentLabel);
//					sibHeight = parentAssigntmentNode.getHeight();
//				}
				srcNodeGrandparent = srcNodeParent.getParent();
				double L = srcNodeGrandparent.getHeight()
						- Math.max(sibHeight, srcNode.getHeight());
				logHR += Math.log(1.0 / L);
			}

			// Detach from current position

			srcNodeParent.removeChild(srcNodeSibling);

			if (srcNodeParent.isRoot()) {
				srcNodeSibling.setParent(null);
			} else {
				srcNodeGrandparent.removeChild(srcNodeParent);
				srcNodeGrandparent.addChild(srcNodeSibling);
			}

			srcNodeParent.setParent(null);
			
			// Define the remaining subtree root
			Node remainingSubtreeRoot;
			if (srcNodeSibling.isRoot())
				remainingSubtreeRoot = srcNodeSibling;
			else
				remainingSubtreeRoot = geneTree.getRoot();

			// Get all nodes, which:
			// 1. are assigned to a transmission node, which is in a subtree rooted at
			// chosenRecipient;
			// 2. are lower than chosenRecipient parent height (transmission event);
			// 3. their parent nodes are higher than chosenRecipient parent height
			// (transmission event);
			fitToAttach = getEdgesAtTransmission(remainingSubtreeRoot, chosenRecipient, recipientGroupNr);
			if (fitToAttach.size() == 0)
				return Double.NEGATIVE_INFINITY;

			// select node at which to attach
			Node newAttachNode;
			newAttachNode = fitToAttach.get(Randomizer.nextInt(fitToAttach.size()));
			
//			int trNodeNr = geneNodeAssignment.get(srcNode.getNr());
//			List<Integer> possibleAssignments = new ArrayList<>();
//			possibleAssignments.add(trNodeNr);
//			possibleAssignments.addAll(getAllParentNrs(intervals.transmissionTreeInput.get().getNode(trNodeNr)));
//			possibleAssignments.addAll(getChildNrs(intervals.transmissionTreeInput.get().getNode(trNodeNr)));

			// store the nodes to choose from for backward move attachment
			fitToAttachAfter = getFitToUniformAttach(remainingSubtreeRoot, srcNode, possibleAssignments);
			if (fitToAttach.size() == 0)
				return Double.NEGATIVE_INFINITY;

			// attach at new height

			srcNodeParent.setHeight(chosenRecipient.getParent().getHeight());

			if (newAttachNode.isRoot()) {
				srcNodeParent.addChild(newAttachNode);
			} else {
				Node oldParent = newAttachNode.getParent();
				oldParent.removeChild(newAttachNode);
				oldParent.addChild(srcNodeParent);
				srcNodeParent.addChild(newAttachNode);
			}


			// Ensure correct root if this has been modified:
			if (srcNodeSibling.isRoot())
				geneTree.setRoot(srcNodeSibling);
			else if (srcNodeParent.isRoot())
				geneTree.setRoot(srcNodeParent);


			if (!Tools.fillAssignmentAndCheck(intervals.transmissionTreeInput.get(), geneTree.getRoot(),
					geneNodeAssignment, trNodeOccupancy))
				return Double.NEGATIVE_INFINITY; // not compatible

			// store the nodes that can be moved in backward move
			List<Node> trueNodesAfter = Pitchforks.getTrueNodes(geneTree);
			fitToMoveAfter = getFitToMoveAtTransmission(trueNodesAfter, chosenRecipient);
			if (fitToMoveAfter.size() == 0)
				return Double.NEGATIVE_INFINITY;

		} else {
			// detach from transmission
//			System.out.println("detach at: " + chosenRecipient.getParent().getHeight());


			fitToMove = getFitToMoveAtTransmission(trueNodes, chosenRecipient);
			if (fitToMove.size() == 0)
				return Double.NEGATIVE_INFINITY;
			Node nodeToMove = fitToMove.get(Randomizer.nextInt(fitToMove.size()));
			Node nodeToMoveParent = nodeToMove.getParent();
			Node sibling = getOtherChild(nodeToMoveParent, nodeToMove);

			// detach from current position

			nodeToMoveParent.removeChild(sibling);

			Node grandParent = null;
			if (nodeToMoveParent.isRoot()) {
				sibling.setParent(null);
			} else {
				grandParent = nodeToMoveParent.getParent();
				grandParent.removeChild(nodeToMoveParent);
				grandParent.addChild(sibling);
			}

			nodeToMoveParent.setParent(null);

			Node remainingSubtreeRoot;
			if (sibling.isRoot())
				remainingSubtreeRoot = sibling;
			else
				remainingSubtreeRoot = geneTree.getRoot();

//			int trNodeNr = geneNodeAssignment.get(nodeToMove.getNr());
//			List<Integer> possibleAssignments = new ArrayList<>();
//			possibleAssignments.add(trNodeNr);
//			possibleAssignments.addAll(getAllParentNrs(intervals.transmissionTreeInput.get().getNode(trNodeNr)));
//			possibleAssignments.addAll(getChildNrs(intervals.transmissionTreeInput.get().getNode(trNodeNr)));

			fitToAttach = getFitToUniformAttach(remainingSubtreeRoot, nodeToMove, possibleAssignments);
			if (fitToAttach.size() == 0)
				return Double.NEGATIVE_INFINITY;
			Node newAttachNode = fitToAttach.get(Randomizer.nextInt(fitToAttach.size()));

			double newHeight;

			if (newAttachNode.isRoot()) {
				double minAttachheight = findMinHeight(newAttachNode, possibleAssignments);
				double offset = Math.max(minAttachheight, nodeToMove.getHeight());
//				if (possibleLocations.contains(geneNodeAssignment.get(newAttachNode.getNr())))
//					offset = Math.max(newAttachNode.getHeight(), nodeToMove.getHeight());
//				else {
//					int parentAssignmentLabel = geneNodeAssignment.get(newAttachNode.getParent().getNr());
//					Node parentAssigntmentNode = intervals.transmissionTreeInput.get().getNode(parentAssignmentLabel);
//					offset = Math.max(parentAssigntmentNode.getHeight(), nodeToMove.getHeight());
//				}
//				double offset = Math.max(nodeToMove.getHeight(), newAttachNode.getHeight());
				double expRate = 1.0 / (rootAttachLambda * offset);
				newHeight = offset + Randomizer.nextExponential(expRate);

				logHR -= -expRate * (newHeight - offset)
						+ Math.log(expRate);
			} else {
				double sibHeight = findMinHeight(newAttachNode, possibleAssignments);
//				if (possibleLocations.contains(geneNodeAssignment.get(newAttachNode.getNr())))
//					sibHeight = newAttachNode.getHeight();
//				else {
//					int parentAssignmentLabel = geneNodeAssignment.get(newAttachNode.getParent().getNr());
//					Node parentAssigntmentNode = intervals.transmissionTreeInput.get().getNode(parentAssignmentLabel);
//					sibHeight = parentAssigntmentNode.getHeight();
//				}

				double L = newAttachNode.getParent().getHeight() -
						Math.max(nodeToMove.getHeight(), sibHeight);
				newHeight = Randomizer.nextDouble() * L +
						Math.max(nodeToMove.getHeight(), sibHeight);
				logHR -= Math.log(1.0 / L);
			}

			fitToAttachAfter = getEdgesAtTransmission(remainingSubtreeRoot, chosenRecipient, recipientGroupNr);

			if (fitToAttachAfter.size() == 0)
				return Double.NEGATIVE_INFINITY;

			// attach at new height

			nodeToMoveParent.setHeight(newHeight);

			if (newAttachNode.isRoot()) {
				nodeToMoveParent.addChild(newAttachNode);
			} else {
				Node oldParent = newAttachNode.getParent();
				oldParent.removeChild(newAttachNode);
				oldParent.addChild(nodeToMoveParent);
				nodeToMoveParent.addChild(newAttachNode);
			}


			// Ensure correct root if set if this has been modified:
			if (sibling.isRoot())
				geneTree.setRoot(sibling);
			else if (nodeToMoveParent.isRoot())
				geneTree.setRoot(nodeToMoveParent);


			if (!Tools.fillAssignmentAndCheck(intervals.transmissionTreeInput.get(), geneTree.getRoot(),
					geneNodeAssignment, trNodeOccupancy))
				return Double.NEGATIVE_INFINITY; // not compatible

			List<Node> trueNodesAfter = Pitchforks.getTrueNodes(geneTree);
//			List<Node> trueNodesWithParentsNotAtTransmissionAfter = Tools.getGeneNodesWithParentNotAtTransmission(
//					trueNodesAfter,
//					trHeights);
			fitToMoveAfter = getNodesToMoveToTransmission(trueNodesAfter,
					chosenRecipient,
					recipientGroupNr);
			if (fitToMoveAfter.size() == 0)
				return Double.NEGATIVE_INFINITY;

		}

		logHR -= Math.log(1.0 / fitToMove.size());
		logHR -= Math.log(1.0 / fitToAttach.size());

		logHR += Math.log(1.0 / fitToMoveAfter.size());
		logHR += Math.log(1.0 / fitToAttachAfter.size());

//		if (fitToMove.size() > 2 || fitToAttach.size() > 2 || fitToMoveAfter.size() > 2 || fitToAttachAfter.size() > 2)
//			System.out.println();

		return logHR;
    }


	private void getRecipients(Node subRoot, List<Node> recipients) {
		if (subRoot.isLeaf())
			return;
		Node leftChild = subRoot.getChild(0);
		getRecipients(leftChild, recipients);
		if (subRoot.getChildCount() > 1 && !subRoot.isFake()) {
			Node rightChild = subRoot.getChild(1);
			recipients.add(rightChild);
			getRecipients(rightChild, recipients);
		}
		return;
    }

	/**
	 * @param trueNodes          list of all logical nodes in a tree
	 * @param trNodeParentHeight height or transmission tree recipient node parent
	 * @return list of true, non-root nodes, which parents are NOT part of polytomy
	 *         or merger and are not at transmission height
	 */
	private List<Node> getNodesToMoveToTransmission(List<Node> trueNodes,
			Node trNode, List<Integer> recipientGroupNr) {
		List<Node> fitToMove = new ArrayList<>();
		List<Node> possibleSubRoots = new ArrayList<>();
		List<Node> allValidNodes = new ArrayList<>();
		Node subRoot = null;

		for (Node n : trueNodes) {
//			if (!n.isRoot() && !Pitchforks.isPolytomy(n) && !Tools.isMultiMerger(trueNodes, n)
			if (recipientGroupNr.contains(geneNodeAssignment[n.getNr()])) {
				if (!recipientGroupNr.contains(geneNodeAssignment[n.getParent().getNr()]))
					possibleSubRoots.add(n);
			if (!n.isRoot()
					&& !Pitchforks.isPolytomy(n.getParent()) && !Tools.isMultiMerger(trueNodes, n.getParent())
						&& !Arrays.stream(trHeights).anyMatch(x -> x == n.getParent().getHeight())// !=
																									// trNode.getParent().getHeight()
						&& n.getHeight() < trNode.getParent().getHeight()) {
				fitToMove.add(n);
			}
		}
		}
		if (possibleSubRoots.size() == 1 && fitToMove.contains(possibleSubRoots.get(0)))
			fitToMove.remove(possibleSubRoots.get(0));


		return fitToMove;
	}

	private List<Node> getFitToMoveAtTransmission(List<Node> trueNodes, Node recipient) {
		List<Node> fitToMove = new ArrayList<>();

		for (Node n : trueNodes) {
//			if (!n.isRoot())
//				if (geneNodeAssignment.get(n.getParent().getNr()) != recipient.getNr()
//						&& n.getParent().getHeight() == recipient.getParent().getHeight())
//					System.out.println();
			if (!n.isRoot()
					&& geneNodeAssignment[n.getParent().getNr()] == recipient.getNr()
					&& n.getParent().getHeight() == recipient.getParent().getHeight())
//					&& Pitchforks.isLogicalNode(n.getParent()))
				fitToMove.add(n);
		}

		return fitToMove;
	}

	/**
	 * Get all nodes, which: 1. are lower than chosenRecipient parent height
	 * (transmission event); 2. their parent nodes are higher than chosenRecipient
	 * parent height (transmission event);
	 * 
	 * @param subtreeRoot
	 * @param recipient
	 * @param recipientGroupNr
	 * @return
	 */
	private List<Node> getEdgesAtTransmission(Node subtreeRoot, Node recipient, List<Integer> recipientGroupNr) {
		List<Node> edgesAtTransmission = new ArrayList<>();

		if (Pitchforks.isLogicalNode(subtreeRoot)
				&& recipientGroupNr.contains(geneNodeAssignment[subtreeRoot.getNr()])
				&& (subtreeRoot.isRoot() ||
						subtreeRoot.getParent().getHeight() > recipient.getParent().getHeight()) // TODO equality here
																									// is an experiment.
																									// CHECK!!
				&&
				subtreeRoot.getHeight() <= recipient.getParent().getHeight()) {
//			if (subtreeRoot.getHeight() == recipient.getParent().getHeight())
//				System.out.println();
			edgesAtTransmission.add(subtreeRoot);
		}

		if (subtreeRoot.getHeight() > recipient.getParent().getHeight()) {
			for (Node child : subtreeRoot.getChildren())
				edgesAtTransmission.addAll(getEdgesAtTransmission(child, recipient, recipientGroupNr));
		}


		return edgesAtTransmission;
	}



	private List<Node> getFitToUniformAttach(Node subtreeRoot, Node nodeToMove, List<Integer> possibleNodeAssignments) {
		List<Node> fitToAttach = new ArrayList<>();

		if (Pitchforks.isLogicalNode(subtreeRoot)
				&& (subtreeRoot.isRoot() || (nodeToMove.getHeight() < subtreeRoot.getParent().getHeight()
						&& possibleNodeAssignments.contains(geneNodeAssignment[subtreeRoot.getParent().getNr()]))))
			fitToAttach.add(subtreeRoot);
		for (Node child : subtreeRoot.getChildren())
				fitToAttach.addAll(getFitToUniformAttach(child, nodeToMove, possibleNodeAssignments));


		return fitToAttach;
	}

	private List<Integer> getAllParentNrs(Node trNode) {
		List<Integer> tmp = new ArrayList<>();
		Node trNodeTmp = trNode;
		while (!trNodeTmp.isRoot()) {
			tmp.add(trNodeTmp.getParent().getNr());
			trNodeTmp = trNodeTmp.getParent();
		}

		return tmp;
	}

	private List<Integer> getChildNrs(Node trNode) {
		List<Integer> tmp = new ArrayList<>();

		if (trNode.isLeaf())
			return tmp;
		for (Node child : trNode.getChildren()) {
			tmp.add(child.getNr());
			tmp.addAll(getChildNrs(child));
		}

		return tmp;
	}

	private double findMinHeight(Node geneNode, List<Integer> possibleNodeAssignments) {
		int parentAssignmentLabel = geneNodeAssignment[geneNode.getNr()];
		if (possibleNodeAssignments.contains(parentAssignmentLabel))
			return geneNode.getHeight();
		Node trNode = intervals.transmissionTreeInput.get().getNode(parentAssignmentLabel);
		while (true) {
			if (possibleNodeAssignments.contains(trNode.getNr()))
				return trNode.getHeight();
			trNode = trNode.getParent();

		}
	}

	/**
	 * @param trNodeNr
	 * @return a single gene tree node, which is embedded in transmission tree
	 *         branch and has max height. It does not mean that it is the only such
	 *         node.
	 */
//	private Node getHeighestGeneNodeInTransmissionBranch(int trNodeNr) {
//		double maxHeight = Double.NEGATIVE_INFINITY;
//		Node maxHeightNode = null;
//
//
//		List<Integer> assignedNodeNrs = inverseGeneNodeAssignment.get(trNodeNr);
//
//		if (assignedNodeNrs != null) {
//			for (int geneNodeNr : assignedNodeNrs) {
//				if (Pitchforks.isLogicalNode(geneTree.getNode(geneNodeNr))
//						&& geneTree.getNode(geneNodeNr).getHeight() > maxHeight) {
//					maxHeightNode = geneTree.getNode(geneNodeNr);
//					maxHeight = maxHeightNode.getHeight();
//				}
//			}
//		}
//
//		return maxHeightNode;
//	}


}
