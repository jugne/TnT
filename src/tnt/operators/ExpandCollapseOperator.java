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
import tnt.distribution.GeneTreeIntervals;
import tnt.util.Tools;

// TODO: Auto-generated Javadoc
/**
 * The Class ExpandCollapseOperator.
 * 
 * 
 */
public class ExpandCollapseOperator extends TreeOperator {

	/** The root attach lambda input. */
    public Input<Double> rootAttachLambdaInput = new Input<>("rootAttachLambda",
            "Mean of exponential (relative to tree height) from which " +
                    "expanded node height is drawn if expanded from a " +
                    "root polytomy.",
            0.1);

	/** The gene tree intervals input. */
	public Input<GeneTreeIntervals> geneTreeIntervalsInput = new Input<>("geneTreeIntervals",
			"intervals for a gene tree", Validate.REQUIRED);

	/** The gene tree. */
    Tree tree;

	/** The lambda for root attachment. */
    double lambda;

	/**
	 * The gene node assignment (gene tree node Nr -> transmission tree node Nr).
	 */
	Integer[] geneNodeAssignment;

	/** The gene tree intervals. */
	GeneTreeIntervals intervals;

	/**
	 * The transmission node heights on the transmission tree (bifurcating node
	 * heights).
	 */
	List<Double> trHeights;

	/**
	 * Initialize the operator
	 */
    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        lambda = rootAttachLambdaInput.get();
		intervals = geneTreeIntervalsInput.get();
    }

	/**
	 * Generate a proposal by this operator.
	 *
	 * @return log of Hastings ratio
	 */
    @Override
    public double proposal() {

        double logHR = 0.0;
		geneNodeAssignment = intervals.getGeneTreeNodeAssignment();
		trHeights = Tools.getTransmissionHeights(intervals.transmissionTreeInput.get());
		boolean which = Randomizer.nextBoolean();

		int p_bf = getPolytomyCount(tree);
		int m_bf = updateMultiMergeCount(tree);
		final Tree treeaf = tree.copy();

		if (which) {
            // Collapse

            List<Node> collapsableEdges = getCollapsableEdges(tree);

            if (collapsableEdges.isEmpty())
                return Double.NEGATIVE_INFINITY;

            Node edgeToCollapse = collapsableEdges.get(Randomizer.nextInt(collapsableEdges.size()));
            logHR -= Math.log(1.0/collapsableEdges.size());

            Node edgeParent = edgeToCollapse.getParent();
            Node sister = getOtherChild(edgeParent, edgeToCollapse);

            if (edgeParent.isRoot()) {
                double expRate = 1.0/(lambda*sister.getHeight());
                logHR += -expRate*(edgeParent.getHeight() - sister.getHeight()) + Math.log(expRate);
            } else {
                double L = edgeParent.getParent().getHeight() - sister.getHeight();
                logHR += Math.log(1.0/L);
            }

            edgeParent.setHeight(sister.getHeight());

            logHR += Math.log(1.0/getExpandableEdges(tree).size());

        } else {
            // Expand

            List<Node> expandableEdges = getExpandableEdges(tree);

            if (expandableEdges.isEmpty())
                return Double.NEGATIVE_INFINITY;

            Node edgeToExpand = expandableEdges.get(Randomizer.nextInt(expandableEdges.size()));
            logHR -= Math.log(1.0/expandableEdges.size());

            Node logicalParent = Pitchforks.getLogicalParent(edgeToExpand);
            assert logicalParent != null;

            double newHeight;
            if (logicalParent.isRoot()) {
                double expRate = 1.0/(lambda*logicalParent.getHeight());
                newHeight = logicalParent.getHeight() + Randomizer.nextExponential(expRate);
                logHR -= -expRate*(newHeight - logicalParent.getHeight()) + Math.log(expRate);
            } else {
                double L = logicalParent.getParent().getHeight() - logicalParent.getHeight();
                newHeight = logicalParent.getHeight() + Randomizer.nextDouble()*L;
                logHR -= Math.log(1.0/L);
            }

            // Disconnect edge

            Node nodeToMove = edgeToExpand.getParent();
            Node sisterNode = getOtherChild(nodeToMove, edgeToExpand);
            Node nodeToMoveParent = nodeToMove.getParent();

            nodeToMove.removeChild(sisterNode);
            if (nodeToMoveParent != null) {
                nodeToMoveParent.removeChild(nodeToMove);
                nodeToMoveParent.addChild(sisterNode);
            } else {
                sisterNode.setParent(null);
            }
            nodeToMove.setParent(null);


			if (nodeToMove.getNr() == logicalParent.getNr())
                logicalParent = sisterNode;

            // Attach edge

            if (sisterNode.isRoot()) {
                nodeToMove.addChild(sisterNode);
                tree.setRoot(nodeToMove);
             } else {
                if (logicalParent.isRoot()) {
                    nodeToMove.addChild(logicalParent);
                    tree.setRoot(nodeToMove);
                } else {
                    Node logicalParentParent = logicalParent.getParent();
                    logicalParentParent.removeChild(logicalParent);
                    logicalParentParent.addChild(nodeToMove);
                    nodeToMove.addChild(logicalParent);
                }
            }

			// Set new node height
			nodeToMove.setHeight(newHeight);

            // Complete HR calculation
			List<Node> colAfter = getCollapsableEdges(tree);

			logHR += Math.log(1.0 / colAfter.size());
        }


		int p_af = getPolytomyCount(tree);
		int m_af = updateMultiMergeCount(tree);

		if (Math.abs(p_bf - p_af) != 1 && Math.abs(p_bf - p_af) != 0) {
			System.out.println();
		}
		if (m_bf != m_af) {
			System.out.println();
		}


        return logHR;
    }

	private int updateMultiMergeCount(Tree tree) {
		// Zero entries
		int nMultiMerge = 0;

		// Compute histogram
		List<Node> trueNodes = Pitchforks.getTrueInternalNodes(tree);
		HashMap<Double, Integer> mergerMap = new HashMap<Double, Integer>();

		for (Node node : trueNodes) {
			if (!mergerMap.keySet().contains(node.getHeight())) {
				mergerMap.put(node.getHeight(), 0);
			} else
				mergerMap.put(node.getHeight(), mergerMap.get(node.getHeight()) + 1);
		}

		for (Double key : mergerMap.keySet()) {
			if (mergerMap.get(key) > 0)
				nMultiMerge += 1;
		}
		return nMultiMerge;
	}

	private int getPolytomyCount(Tree tree) {
		int count = 0;

		List<Node> trueNodes = new ArrayList<>();
		for (Node node : tree.getNodesAsArray())
			if (node.isRoot() || node.getParent().getHeight() > node.getHeight())
				trueNodes.add(node);

		for (Node node : trueNodes) {
			if (!node.isLeaf() && (node.getChildren().get(0).getHeight() == node.getHeight()
					|| node.getChildren().get(1).getHeight() == node.getHeight()))
				count += 1;
		}

		return count;
	}

	/**
	 * Gets the collapsable edges.
	 *
	 * @param tree the gene tree
	 * @return the collapsable edges of this gene tree.
	 * 
	 *         Edge is "collapsable" if: 1. parent node is NOT at transmission event
	 *         (technically difficult to manage such move) 2.
	 */
    private List<Node> getCollapsableEdges(Tree tree) {

        List<Node> trueNodes = Pitchforks.getTrueNodes(tree);
        List<Node> collapsableEdges = new ArrayList<>();
		List<Node> noTrParentNodes = Tools.getGeneNodesWithParentNotAtTransmission(trueNodes, trHeights);


		for (Node node : noTrParentNodes) {
			int trNodeNr = geneNodeAssignment[node.getNr()];
			List<Integer> possibleAssignments = new ArrayList<>();
			possibleAssignments.add(trNodeNr);
			possibleAssignments.addAll(Tools.getAllParentNrs(intervals.transmissionTreeInput.get().getNode(trNodeNr)));
			possibleAssignments.addAll(Tools.getChildNrs(intervals.transmissionTreeInput.get().getNode(trNodeNr)));

			if (node.isRoot() || Pitchforks.isPolytomy(node.getParent())
					|| Tools.isMultiMerger(trueNodes, node.getParent()))
                    continue;

            Node sister = getOtherChild(node.getParent(), node);
			int trNodeSisterNr = geneNodeAssignment[sister.getNr()];
			if (Tools.greaterOrEqualHeightNode(node, sister) || sister.isLeaf()
					|| !possibleAssignments.contains(trNodeSisterNr))
                continue;

            collapsableEdges.add(node);
        }

        return collapsableEdges;
    }

	/**
	 * Gets the expandable edges.
	 *
	 * @param tree the tree
	 * @return the expandable edges
	 */
    private List<Node> getExpandableEdges(Tree tree) {

        List<Node> trueNodes = Pitchforks.getTrueNodes(tree);
        List<Node> expandableEdges = new ArrayList<>();

        for (Node node : trueNodes) {
			if (!node.isRoot() && Pitchforks.isPolytomy(node.getParent()))
                expandableEdges.add(node);
        }

        return expandableEdges;
    }
}
