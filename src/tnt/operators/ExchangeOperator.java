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
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import pitchfork.Pitchforks;
import tnt.distribution.GeneTreeIntervals;
import tnt.util.Tools;

@Description("Exchange operator compatible with pitchfork trees.")
public class ExchangeOperator extends TreeOperator {

    public Input<Boolean> isNarrowInput = new Input<>(
            "isNarrow",
            "Whether narrow exchange is used.",
            true);

	/** The gene tree intervals input. */
	public Input<GeneTreeIntervals> geneTreeIntervalsInput = new Input<>("geneTreeIntervals",
			"intervals for a gene tree", Validate.REQUIRED);

    boolean isNarrow;
	/** The gene tree. */
    Tree tree;

	/**
	 * The gene node assignment (gene tree node Nr -> transmission tree node Nr).
	 */
	Integer[] geneNodeAssignment;

	/** The gene tree intervals. */
	GeneTreeIntervals intervals;


    @Override
    public void initAndValidate() {
        isNarrow = isNarrowInput.get();
        tree = treeInput.get();
		intervals = geneTreeIntervalsInput.get();
    }

    @Override
    public double proposal() {

		geneNodeAssignment = intervals.getGeneTreeNodeAssignment();

        if (isNarrow) {
            List<Node> trueNodes = Pitchforks.getTrueNodes(tree);

            if (trueNodes.size() - tree.getLeafNodeCount() <= 1)
                return Double.NEGATIVE_INFINITY;

            Node srcNode, srcNodeParent, destNode, destNodeParent;
            Node srcNodeLogicalParent, srcNodeLogicalGrandparent = null;

            do {
                srcNode = trueNodes.get(Randomizer.nextInt(trueNodes.size()));
                srcNodeLogicalParent = Pitchforks.getLogicalParent(srcNode);

                if (srcNodeLogicalParent != null)
                    srcNodeLogicalGrandparent = Pitchforks.getLogicalParent(srcNodeLogicalParent);

            } while (srcNodeLogicalParent == null || srcNodeLogicalGrandparent == null);
            srcNodeParent = srcNode.getParent();

			List<Node> possibleDestNodes = getPossibleDestNodes(srcNodeLogicalGrandparent, srcNodeLogicalParent,
					srcNode);
			if (possibleDestNodes.isEmpty())
				return Double.NEGATIVE_INFINITY;

            do {
                destNode = possibleDestNodes.get(Randomizer.nextInt(possibleDestNodes.size()));
            } while (destNode == srcNodeLogicalParent);
            destNodeParent = destNode.getParent();

            // Reject if substitution would result in negative edge length:
			if (destNode.getHeight() >= srcNodeParent.getHeight()
					|| srcNode.getHeight() >= destNodeParent.getHeight())
                return Double.NEGATIVE_INFINITY;

            srcNodeParent.removeChild(srcNode);
            destNodeParent.removeChild(destNode);
            srcNodeParent.addChild(destNode);
            destNodeParent.addChild(srcNode);

        } else {
            throw new UnsupportedOperationException("Wide exchange for " +
                    "pitchfork trees is not yet supported.");
        }


        return 0;
    }

	private List<Node> getPossibleDestNodes(Node srcNodeLogicalGrandparent, Node srcNodeLogicalParent, Node srcNode) {
		List<Node> possibleDestNodes = Pitchforks.getLogicalChildren(srcNodeLogicalGrandparent);
		List<Integer> possibleAssignments = new ArrayList<>();
		int trNodeNr = geneNodeAssignment[srcNode.getNr()];
		possibleAssignments.add(trNodeNr);
		possibleAssignments.addAll(Tools.getAllParentNrs(intervals.transmissionTreeInput.get().getNode(trNodeNr)));
		possibleAssignments.addAll(Tools.getChildNrs(intervals.transmissionTreeInput.get().getNode(trNodeNr)));
		
		List<Node> remove = new ArrayList<Node>();
		for (Node n : possibleDestNodes) {
			if (!possibleAssignments.contains(geneNodeAssignment[n.getNr()])
					|| Tools.equalHeightWithPrecisionNode(n, srcNodeLogicalParent))
					remove.add(n);
		}
		
		possibleDestNodes.removeAll(remove);
		
		return possibleDestNodes;

	}
}
