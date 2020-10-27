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
import java.util.HashMap;
import java.util.List;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import pitchfork.Pitchforks;
import tnt.distribution.GeneTreeEvent;
import tnt.distribution.GeneTreeIntervals;
import tnt.util.Tools;

public class TransmissionAttach4 extends TreeOperator {

	public Input<GeneTreeIntervals> geneTreeIntervalsInput = new Input<>("geneTreeIntervals",
			"intervals for a gene tree", Validate.REQUIRED);

	public Input<Double> rootAttachLambdaInput = new Input<>(
			"rootAttachLambda",
			"Mean of exponential distribution (relative to tree height)" +
					"used to position attachments above the root.",
			2.0);

    Tree tree;
	GeneTreeIntervals intervals;
	HashMap<Integer, List<GeneTreeEvent>> eventsPerTransmissionTreeNode;
	HashMap<Integer, Integer> geneNodeAssignment;
	Double rootAttachLambda;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
		intervals = geneTreeIntervalsInput.get();
		rootAttachLambda = rootAttachLambdaInput.get();

    }

    @Override
    public double proposal() {

        double logHR = 0.0;
        eventsPerTransmissionTreeNode = intervals.getGeneTreeEventList();
		geneNodeAssignment = intervals.getGeneTreeNodeAssignment();
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
		for (Node n : chosenRecipient.getAllChildNodesAndSelf()) {
			recipientGroupNr.add(n.getNr());
		}

		List<Node> trueNodes = Pitchforks.getTrueNodes(tree);
		List<Node> fitToMove = new ArrayList<>();
		List<Node> fitToMoveAfter = new ArrayList<>();
		List<Node> fitToAttach = new ArrayList<>();
		List<Node> fitToAttachAfter = new ArrayList<>();

		if (Randomizer.nextBoolean()) {
			// attach a node at transmission
			
			// get all nodes, that are not at transmission height and can be attached there
			fitToMove = getNodesToMoveToTransmission(trueNodes, chosenRecipient.getParent().getHeight());
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
				double offset = Math.max(srcNodeSibling.getHeight(), srcNode.getHeight());
				double expRate = 1.0 / (rootAttachLambda * offset);
				logHR += -expRate * (srcNodeParent.getHeight() - offset) + Math.log(expRate);
			} else {
				srcNodeGrandparent = srcNodeParent.getParent();
				double L = srcNodeGrandparent.getHeight()
						- Math.max(srcNodeSibling.getHeight(), srcNode.getHeight());
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
				remainingSubtreeRoot = tree.getRoot();

			// Get all nodes, which:
			// 1. are assigned to a transmission node, which is in a subtree rooted at
			// chosenRecipient;
			// 2. are lower than chosenRecipient parent height (transmission event);
			// 3. their parent nodes are heigher than chosenRecipient parent height
			// (transmission event);
			fitToAttach = getEdgesAtTransmission(remainingSubtreeRoot, chosenRecipient, recipientGroupNr);
			if (fitToAttach.size() == 0)
				return Double.NEGATIVE_INFINITY;

			// select node at which to attach
			Node newAttachNode;
			newAttachNode = fitToAttach.get(Randomizer.nextInt(fitToAttach.size()));

			// store the nodes to choose from for backward move attachment
			fitToAttachAfter = getFitToUniformAttach(remainingSubtreeRoot, srcNode);
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

			// Ensure correct root if set if this has been modified:
			if (srcNodeSibling.isRoot())
				tree.setRoot(srcNodeSibling);
			else if (srcNodeParent.isRoot())
				tree.setRoot(srcNodeParent);

			// store the nodes that can be moved in backward move
			List<Node> trueNodesAfter = Pitchforks.getTrueNodes(tree);
			fitToMoveAfter = getFitToMoveAtTransmission(trueNodesAfter, chosenRecipient);
			if (fitToMoveAfter.size() == 0)
				return Double.NEGATIVE_INFINITY;


		} else {
			// detach from transmission
			
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
				remainingSubtreeRoot = tree.getRoot();

			fitToAttach = getFitToUniformAttach(remainingSubtreeRoot, nodeToMove);
//			fitToAttach.remove(nodeToMoveParent);
//			List<Node> nodesToRemove = Pitchforks.getLogicalChildren(nodeToMove);
//			fitToAttach.removeAll(nodesToRemove);
			if (fitToAttach.size() == 0)
				return Double.NEGATIVE_INFINITY;
			Node newAttachNode = fitToAttach.get(Randomizer.nextInt(fitToAttach.size()));

			double newHeight;

			if (newAttachNode.isRoot()) {
				double offset = Math.max(nodeToMove.getHeight(), newAttachNode.getHeight());
				double expRate = 1.0 / (rootAttachLambda * offset);
				newHeight = offset + Randomizer.nextExponential(expRate);

				logHR -= -expRate * (newHeight - offset)
						+ Math.log(expRate);
			} else {
				double L = newAttachNode.getParent().getHeight() -
						Math.max(nodeToMove.getHeight(), newAttachNode.getHeight());
				newHeight = Randomizer.nextDouble() * L +
						Math.max(nodeToMove.getHeight(), newAttachNode.getHeight());

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
				tree.setRoot(sibling);
			else if (nodeToMoveParent.isRoot())
				tree.setRoot(nodeToMoveParent);

			List<Node> trueNodesAfter = Pitchforks.getTrueNodes(tree);
			fitToMoveAfter = getNodesToMoveToTransmission(trueNodesAfter, chosenRecipient.getParent().getHeight());
			if (fitToMoveAfter.size() == 0)
				return Double.NEGATIVE_INFINITY;


		}

		logHR -= Math.log(1.0 / fitToMove.size());
		logHR -= Math.log(1.0 / fitToAttach.size());

		logHR += Math.log(1.0 / fitToMoveAfter.size());
		logHR += Math.log(1.0 / fitToAttachAfter.size());

		if (tree.getLeafNodeCount() < 16)
			System.exit(-1);

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
			double trNodeParentHeight) {
		List<Node> fitToMove = new ArrayList<>();

		for (Node n : trueNodes) {
//			if (!n.isRoot() && !Pitchforks.isPolytomy(n) && !Tools.isMultiMerger(trueNodes, n)
			if (!n.isRoot() && !Pitchforks.isPolytomy(n.getParent()) && !Tools.isMultiMerger(trueNodes, n.getParent())
					&& n.getParent().getHeight() != trNodeParentHeight && n.getHeight() != trNodeParentHeight) {
				fitToMove.add(n);
			}
		}
		return fitToMove;
	}

	private List<Node> getFitToMoveAtTransmission(List<Node> trueNodes, Node recipient) {
		List<Node> fitToMove = new ArrayList<>();

		for (Node n : trueNodes) {
			if (!n.isRoot() && geneNodeAssignment.get(n.getParent().getNr()) == recipient.getNr()
					&& n.getParent().getHeight() == recipient.getParent().getHeight())
				fitToMove.add(n);
		}

		return fitToMove;
	}

	private List<Node> getEdgesAtTransmission(Node subtreeRoot, Node recipient, List<Integer> recipientGroupNr) {
		List<Node> edgesAtTransmission = new ArrayList<>();

		if (Pitchforks.isLogicalNode(subtreeRoot)
				&& recipientGroupNr.contains(geneNodeAssignment.get(subtreeRoot.getNr()))
				&& (subtreeRoot.isRoot() || subtreeRoot.getParent().getHeight() > recipient.getParent().getHeight()))
			edgesAtTransmission.add(subtreeRoot);

		if (subtreeRoot.getHeight() > recipient.getParent().getHeight()) {
			for (Node child : subtreeRoot.getChildren())
				edgesAtTransmission.addAll(getEdgesAtTransmission(child, recipient, recipientGroupNr));
		}


		return edgesAtTransmission;
	}



	private List<Node> getFitToUniformAttach(Node subtreeRoot, Node nodeToMove) {
		List<Node> fitToAttach = new ArrayList<>();

		if (Pitchforks.isLogicalNode(subtreeRoot) && nodeToMove.getHeight() < subtreeRoot.getHeight()) {
			fitToAttach.add(subtreeRoot);
			for (Node child : subtreeRoot.getChildren())
				fitToAttach.addAll(getFitToUniformAttach(child, nodeToMove));
		}

		return fitToAttach;
	}


}
