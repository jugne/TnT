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
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import pitchfork.Pitchforks;
import tnt.distribution.GeneTreeEvent;
import tnt.distribution.GeneTreeIntervals;

@Description("SPR operator for trees with polytomies.")
public class SPROperator extends TreeOperator {

	public Input<Double> rootAttachLambdaInput = new Input<>(
			"rootAttachLambda",
            "Mean of exponential distribution (relative to tree height)" +
                    "used to position attachments above the root.",
            2.0);

    public Input<Double> probCoalAttachInput = new Input<>(
            "probCoalAttach",
            "Probability of attaching to existing coalescent event.",
            0.1);

	public Input<Double> probMultiMergerInput = new Input<>(
			"probMultiMerger",
			"Probability of attaching at an existing coalescent event height.",
			0.1);

	public Input<GeneTreeIntervals> geneTreeIntervalsInput = new Input<>("geneTreeIntervals",
			"intervals for a gene tree", Validate.REQUIRED);

    Tree tree;
	Double rootAttachLambda, probCoalAttach, probMultiMerger;
	GeneTreeIntervals intervals;
	HashMap<Integer, Integer> geneTreeNodeAssignment;

	HashMap<Integer, List<Integer>> logicalGeneNodesPerTransmissionNode;
	HashMap<Integer, List<GeneTreeEvent>> eventsPerTransmissionTreeNode;

    @Override
    public void initAndValidate() {
        rootAttachLambda = rootAttachLambdaInput.get();
        tree = treeInput.get();
        probCoalAttach = probCoalAttachInput.get();
		probMultiMerger = probMultiMergerInput.get();
		intervals = geneTreeIntervalsInput.get();

		if (probCoalAttach + probMultiMerger > 1) {
			System.err.print("Sum of probabilities is larger than one: probCoalAttach=" + probCoalAttach
					+ ", probMultiMerger=" + probMultiMerger);
			System.exit(1);
		}

    }

    @Override
    public double proposal() {

		geneTreeNodeAssignment = intervals.getGeneTreeNodeAssignment();

		logicalGeneNodesPerTransmissionNode = intervals.getLogicalGeneNodesPerTransmissionNode();
		eventsPerTransmissionTreeNode = intervals.getGeneTreeEventList();

        double logHR = 0.0;

        // Get list of nodes below finite-length edges
        List<Node> trueNodes = getTrueNodes(tree);

        // Record number of (true) edges in original tree:
        int nEdges = trueNodes.size() - 1;
		int nOrigMultiMergerCandidates = 0;
		int nMultiMergerCandidates = 0;

        // Select non-root subtree at random

        Node srcNode, srcNodeParent;
        do {
            srcNode = trueNodes.get(Randomizer.nextInt(trueNodes.size()));
            srcNodeParent = srcNode.getParent();
        } while (srcNodeParent == null);

        Node srcNodeSister = getOtherChild(srcNodeParent, srcNode);

		// Get transmission node that this node is embedded into
		Node trNode = intervals.transmissionTreeInput.get().getNode(geneTreeNodeAssignment.get(srcNode.getNr()));
		List<GeneTreeEvent> eventsInTransmissionBranch = eventsPerTransmissionTreeNode.get(trNode.getNr());
		// Record whether the the original attachment was a multiple merger


        // Record whether the the original attachment was a polytomy
        boolean origAttachWasPolytomy = isPolytomy(srcNodeParent);

		// Record whether the the original attachment was a multiple merger
		// Do this only if origAttach was not polytomy
		// This is because of it can be a multiple merger from before, but only formed a
		// polytomy by attaching to existing node with probCoalAttach. That is, this
		// event would have been more recent and reversing would return us to the start
		// state.
		boolean origAttachWasMultiMerger = false;
		if (!origAttachWasPolytomy) {
			for (GeneTreeEvent e : eventsInTransmissionBranch) {
				if (e.multiCoalSize.size() > 1 && e.time == srcNodeParent.getHeight()) {
					origAttachWasMultiMerger = true;
				}

			}
		}

//        // Incorporate probability of existing attachment point into HR
//
//		if (origAttachWasPolytomy) {
//			logHR += Math.log(probCoalAttach);
//		} else if (origAttachWasMultiMerger) {
//			logHR += Math.log(probMultiMerger);
//		} else {
//            if (!srcNodeSister.isLeaf() && srcNodeSister.getHeight() > srcNode.getHeight())
//				logHR += Math.log(1 - probCoalAttach);
//
//            if (srcNodeParent.isRoot()) {
//                double offset = Math.max(srcNodeSister.getHeight(), srcNode.getHeight());
//                double expRate = 1.0/(rootAttachLambda*offset);
//                logHR += -expRate*(srcNodeParent.getHeight() - offset) + Math.log(expRate);
//            } else {
//                double L = srcNodeParent.getParent().getHeight() - Math.max(srcNodeSister.getHeight(), srcNode.getHeight());
//                logHR += Math.log(1.0/L);
//            }
//        }

		boolean parentWasRoot = srcNodeParent.isRoot();
		Node srcNodeGrandparent = null;

        // Disconnect subtree

        srcNodeParent.removeChild(srcNodeSister);

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
        Node attachmentNode = subtreeNodes.get(Randomizer.nextInt(subtreeNodes.size()));
		List<Node> multiMergSisterCandidates = getNodesForMultiMerger(subtreeNodes, attachmentNode,
				srcNode.getHeight());
		Node multiMergeSister = null;
		nMultiMergerCandidates = multiMergSisterCandidates.size();
		if (nMultiMergerCandidates > 0)
			multiMergeSister = multiMergSisterCandidates.get(Randomizer.nextInt(multiMergSisterCandidates.size()));
		if (origAttachWasMultiMerger) {
			List<Node> origMultiMergSisterCandidates = getNodesForMultiMerger(subtreeNodes, srcNodeSister,
					srcNode.getHeight());
			nOrigMultiMergerCandidates = origMultiMergSisterCandidates.size();

		}
		
        // Incorporate probability of existing attachment point into HR

		if (origAttachWasPolytomy) {
			logHR += Math.log(probCoalAttach);
		} else if (origAttachWasMultiMerger) {
			logHR += Math.log(probMultiMerger);
		} else {
            if (!srcNodeSister.isLeaf() && srcNodeSister.getHeight() > srcNode.getHeight()) {
            	if (nOrigMultiMergerCandidates == 0)
            		logHR += Math.log(1 - probCoalAttach);
            	else if (nOrigMultiMergerCandidates > 0)
            		logHR += Math.log(1 - probCoalAttach - probMultiMerger);
			} else if (nOrigMultiMergerCandidates > 0)
				logHR += Math.log(1 - probMultiMerger);
				

			if (parentWasRoot) {
                double offset = Math.max(srcNodeSister.getHeight(), srcNode.getHeight());
                double expRate = 1.0/(rootAttachLambda*offset);
                logHR += -expRate*(srcNodeParent.getHeight() - offset) + Math.log(expRate);
            } else {
				double L = srcNodeGrandparent.getHeight() - Math.max(srcNodeSister.getHeight(), srcNode.getHeight());
                logHR += Math.log(1.0/L);
            }
        }

        // Determine whether polytomy is to be created

        boolean newAttachIsPolytomy;
		boolean newAttachIsMultiMerger;

        if (attachmentNode.isLeaf() || attachmentNode.getHeight() < srcNode.getHeight()) {
            newAttachIsPolytomy = false;
			if (multiMergeSister != null) {// && !origAttachWasPolytomy) { // no reversibility for polytomies currently
				if (Randomizer.nextDouble() < probMultiMerger) {
					newAttachIsMultiMerger = true;
					logHR -= Math.log(probMultiMerger);
				} else {
					newAttachIsMultiMerger = false;
					logHR -= Math.log(1 - probMultiMerger);
				}
			} else {
				newAttachIsMultiMerger = false;
            }
		} else if (multiMergeSister == null) {// || origAttachWasPolytomy) { // no reversibility for polytomies
												// currently
			newAttachIsMultiMerger = false;

			if (Randomizer.nextDouble() < probCoalAttach) {
				newAttachIsPolytomy = true;
				logHR -= Math.log(probCoalAttach);
			} else {
				newAttachIsPolytomy = false;
				logHR -= Math.log(1 - probCoalAttach);
			}
		} else {
			int s = Randomizer.randomChoicePDF(
					new double[] { probCoalAttach, probMultiMerger, 1 - probCoalAttach - probMultiMerger });

			if (s == 0) {
				newAttachIsPolytomy = true;
				newAttachIsMultiMerger = false;
                logHR -= Math.log(probCoalAttach);
			} else if (s == 1) {
				newAttachIsMultiMerger = true;
				newAttachIsPolytomy = false;
				logHR -= Math.log(probMultiMerger);
            } else {
                newAttachIsPolytomy = false;
				newAttachIsMultiMerger = false;
				logHR -= Math.log(1 - probCoalAttach - probMultiMerger);
            }
        }

        // Select new attachment height

        double attachmentHeight;

        if (newAttachIsPolytomy) {
            attachmentHeight = attachmentNode.getHeight();
		} else if (newAttachIsMultiMerger) {
			attachmentHeight = multiMergeSister.getHeight();
        } else {
            if (attachmentNode.isRoot()) {
                double offset = Math.max(srcNode.getHeight(), attachmentNode.getHeight());
                double expRate = 1.0/(rootAttachLambda*offset);
                attachmentHeight = offset + Randomizer.nextExponential(expRate);

                logHR -= -expRate*(attachmentHeight-offset)
                        + Math.log(expRate);
            } else {
                double L = attachmentNode.getParent().getHeight() -
                        Math.max(srcNode.getHeight(), attachmentNode.getHeight());
                attachmentHeight = Randomizer.nextDouble()*L +
                        Math.max(srcNode.getHeight(), attachmentNode.getHeight());

                logHR -= Math.log(1.0/L);
            }
        }

        // Reconnect subtree

        srcNodeParent.setHeight(attachmentHeight);

        if (attachmentNode.isRoot()) {
            srcNodeParent.addChild(attachmentNode);
        } else {
            Node oldParent = attachmentNode.getParent();
            oldParent.removeChild(attachmentNode);
            oldParent.addChild(srcNodeParent);
            srcNodeParent.addChild(attachmentNode);
        }

        // Ensure correct root if set if this has been modified:
        if (srcNodeSister.isRoot())
            tree.setRoot(srcNodeSister);
        else if (srcNodeParent.isRoot())
            tree.setRoot(srcNodeParent);

        // Account for edge selection probability in HR:
        if (origAttachWasPolytomy != newAttachIsPolytomy) {
            if (origAttachWasPolytomy) {
                logHR += Math.log(nEdges/(nEdges+1.0));
            } else {
                logHR += Math.log(nEdges/(nEdges-1.0));
            }
        }

		if (origAttachWasMultiMerger && !origAttachWasPolytomy) {
			logHR += Math.log(1.0 / nOrigMultiMergerCandidates);
		}
		if (newAttachIsMultiMerger && !newAttachIsPolytomy) {
			logHR -= Math.log(1.0 / nMultiMergerCandidates);
		}
		if (logHR == Double.POSITIVE_INFINITY)
			System.out.println();
        return logHR;
    }

    private List<Node> getNodesInSubtree(Node subtreeRoot, double minAge) {
        List<Node> nodeList = new ArrayList<>();

        if (subtreeRoot.isRoot() || subtreeRoot.getParent().getHeight()>subtreeRoot.getHeight())
            nodeList.add(subtreeRoot);

        if (subtreeRoot.getHeight()>minAge) {
            for (Node child : subtreeRoot.getChildren())
                nodeList.addAll(getNodesInSubtree(child, minAge));
        }

        return nodeList;
    }

	private List<Node> getNodesForMultiMerger(List<Node> nodesInSubtree, Node newAttach, double minAge) {
		List<Node> nodeList = new ArrayList<>();
		Set<Double> heights = new HashSet<Double>();
		
		if (newAttach.isRoot())
			return nodeList;
		for (Node n : nodesInSubtree) {
			if (!n.isLeaf() && n != newAttach && n.getHeight() < Pitchforks.getLogicalParent(newAttach).getHeight()
					&& n.getHeight() > newAttach.getHeight() && n.getHeight() > minAge
					&& !heights.contains(n.getHeight())) { // only add possible node heights once
				nodeList.add(n);
				heights.add(n.getHeight());
			}
		}
		return nodeList;
	}
}
