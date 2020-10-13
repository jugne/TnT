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

import static pitchfork.Pitchforks.getGroupAndLogicalChildren;
import static pitchfork.Pitchforks.getTrueInternalNodes;

import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import pitchfork.Pitchforks;

@Description("Uniform node height operator compatible with trees having polytomies.")
public class UniformOperator extends TreeOperator {

    public Input<Boolean> scaleRootInput = new Input<>(
            "scaleRoot",
            "Whether to scale the age of the root node.",
            true);


    public Input<Double> scaleFactorInput = new Input<>(
            "scaleFactor",
            "Tuning parameter for scaling root.",
            0.8);

    Tree tree;

    boolean scaleRoot;
    double scaleFactor;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        scaleRoot = scaleRootInput.get();
        scaleFactor = scaleFactorInput.get();
    }

    @Override
    public double proposal() {

        double logHR = 0.0;

		boolean makeMerger = Randomizer.nextDouble() < 0.9;
		List<Node> trueNodes = getTrueInternalNodes(tree);

		if (makeMerger) {
			logHR -= Math.log(0.9);
			// if one true inner node = everything is one polytomy. cannot operate on that
			// if two inner nodes = one polytomy, plus root. cannot operate on that
			int nTrueInnerNodes = trueNodes.size();
			if (nTrueInnerNodes < 2)
				return Double.NEGATIVE_INFINITY;

			Node srcNode;
			do {
				srcNode = trueNodes.get(Randomizer.nextInt(nTrueInnerNodes));
			} while (srcNode.isRoot());
			logHR -= Math.log(1.0 / (nTrueInnerNodes - 1));

			Node trueParent = Pitchforks.getLogicalParent(srcNode);
			double maxHeight = trueParent.getHeight();
			List<Node> nodesInLogicalGroup = new ArrayList<>();
			List<Node> logicalChildren = new ArrayList<>();
			getGroupAndLogicalChildren(srcNode, nodesInLogicalGroup, logicalChildren);
			double minHeight = logicalChildren.stream().mapToDouble(Node::getHeight).max().getAsDouble();

			List<Double> possibleMergerHeights = new ArrayList<Double>();
			boolean wasSrcNodeInMerger = false;
			for (Node n : trueNodes) {
				if (n.getHeight() == srcNode.getHeight() && n.getNr() != srcNode.getNr())
					wasSrcNodeInMerger = true;
				else if (!n.isRoot() && n.getNr() != srcNode.getNr() && n.getHeight() > minHeight
						&& n.getHeight() < maxHeight) {
					possibleMergerHeights.add(n.getHeight());
				}
			}
			int nNewMergerHeights = possibleMergerHeights.size();
			if (nNewMergerHeights == 0)
				return Double.NEGATIVE_INFINITY; // no nodes to choose from
			double newMergerHeight = possibleMergerHeights.get(Randomizer.nextInt(nNewMergerHeights));
			logHR -= Math.log(1.0 / nNewMergerHeights);

			if (wasSrcNodeInMerger) {
				logHR += Math.log(1.0 / (nTrueInnerNodes - 1));
				logHR += Math.log(1.0 / nNewMergerHeights);
				logHR += Math.log(0.9);
			} else {
				double L = maxHeight - minHeight;
				logHR += Math.log(1.0 / L);
				logHR += Math.log(1.0 / (nTrueInnerNodes - 1)); // srcNode could not have been a root
				logHR += Math.log(0.1);
			}

			srcNode.setHeight(newMergerHeight);
			for (Node node : nodesInLogicalGroup)
				node.setHeight(newMergerHeight);

		} else {
			logHR -= Math.log(0.1);
			int nTrueInnerNodes = trueNodes.size();
			if (nTrueInnerNodes == 1 && !scaleRoot)
				return Double.NEGATIVE_INFINITY;
			Node srcNode;
			do {
				srcNode = trueNodes.get(Randomizer.nextInt(nTrueInnerNodes));
			} while (!scaleRoot && srcNode.isRoot());

			// account for choosing a node
			if (scaleRoot)
				logHR -= Math.log(1.0 / nTrueInnerNodes);
			else
				logHR -= Math.log(1.0 / (nTrueInnerNodes - 1));

			List<Node> nodesInLogicalGroup = new ArrayList<>();
			List<Node> logicalChildren = new ArrayList<>();
			getGroupAndLogicalChildren(srcNode, nodesInLogicalGroup, logicalChildren);
			double minHeight = logicalChildren.stream().mapToDouble(Node::getHeight).max().getAsDouble();

			double newHeight;
			if (srcNode.isRoot()) {

				double minf = Math.min(scaleFactor, 1.0 / scaleFactor);
				double maxf = 1.0 / minf;
				double f = minf + Randomizer.nextDouble() * (maxf - minf);

				newHeight = srcNode.getHeight() * f;

				if (newHeight < minHeight)
					return Double.NEGATIVE_INFINITY;

				logHR = -Math.log(f);
//				logHR += Math.log(1.0 / nTrueInnerNodes);
			} else {
				Node parent = srcNode.getParent();
				double maxHeight = parent.getHeight();
				double L = maxHeight - minHeight;
				newHeight = minHeight + (maxHeight - minHeight) * Randomizer.nextDouble();

				logHR -= Math.log(1.0 / L);

				boolean wasSrcNodeInMerger = false;
				List<Double> possibleMergerHeights = new ArrayList<Double>();
				for (Node n : trueNodes) {
					if (n.getHeight() == srcNode.getHeight() && n.getNr() != srcNode.getNr()) {
						wasSrcNodeInMerger = true;
						possibleMergerHeights.add(n.getHeight());
					}
				}
				if (wasSrcNodeInMerger) {
					int nNewMergerHeights = possibleMergerHeights.size();
					logHR += Math.log(1.0 / nNewMergerHeights);
					logHR += Math.log(1.0 / (nTrueInnerNodes - 1));
					logHR += Math.log(0.9);
				} else {
					logHR += Math.log(1.0 / L);
					logHR += Math.log(0.1);
					if (scaleRoot)
						logHR += Math.log(1.0 / nTrueInnerNodes);
					else
						logHR += Math.log(1.0 / (nTrueInnerNodes - 1));
				}
			}

			srcNode.setHeight(newHeight);
			for (Node node : nodesInLogicalGroup)
				node.setHeight(newHeight);

		}

		return logHR;
    }
}
