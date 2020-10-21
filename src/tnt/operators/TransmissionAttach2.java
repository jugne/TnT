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
import tnt.distribution.GeneTreeEvent.GeneTreeEventType;
import tnt.distribution.GeneTreeIntervals;
import tnt.util.Tools;

public class TransmissionAttach2 extends TreeOperator {

	public Input<GeneTreeIntervals> geneTreeIntervalsInput = new Input<>("geneTreeIntervals",
			"intervals for a gene tree", Validate.REQUIRED);

    Tree tree;
	GeneTreeIntervals intervals;
	HashMap<Integer, List<GeneTreeEvent>> eventsPerTransmissionTreeNode;
	HashMap<Integer, Integer> nodesPerTransmissionTreeNode;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
		intervals = geneTreeIntervalsInput.get();

    }

    @Override
    public double proposal() {

        double logHR = 0.0;
        eventsPerTransmissionTreeNode = intervals.getGeneTreeEventList();
		nodesPerTransmissionTreeNode = intervals.getGeneTreeNodeAssignment();
		// fill the recipient node list on
		List<Node> recipients = new ArrayList<>();
		Node trTreeRoot = intervals.transmissionTreeInput.get().getRoot();
		getRecipients(trTreeRoot, recipients);

		int nRecipients = recipients.size();
		Node chosenRecipient = recipients.get(Randomizer.nextInt(nRecipients));
		List<GeneTreeEvent> eventList = eventsPerTransmissionTreeNode.get(chosenRecipient.getNr());
		if (eventList.size() == 0)
			return Double.NEGATIVE_INFINITY;


        if (Randomizer.nextBoolean()) {
			// Up


			List<Node> fitToMove = getfitNodesUp(eventList, chosenRecipient.getParent().getHeight());

			if (fitToMove.size() == 0)
				return Double.NEGATIVE_INFINITY;
			Node nodeToMove = fitToMove.get(Randomizer.nextInt(fitToMove.size()));
			double heightMin = Math.max(nodeToMove.getChild(0).getHeight(), nodeToMove.getChild(1).getHeight());
			 
			nodeToMove.setHeight(chosenRecipient.getParent().getHeight());

			List<Node> fitToMoveAfter = getFitNodesDown(eventList, chosenRecipient.getParent().getHeight());


			logHR += Math.log(1.0 / fitToMoveAfter.size());
			logHR += Math.log(1.0 / (chosenRecipient.getParent().getHeight() - heightMin));
			logHR -= Math.log(1.0 / fitToMove.size());
			

        } else {
			// Down

//			System.out.println("Down");
			List<Node> fitToMove = getFitNodesDown(eventList, chosenRecipient.getParent().getHeight());
			if (fitToMove.size() == 0)
				return Double.NEGATIVE_INFINITY;
			Node nodeToMove = fitToMove.get(Randomizer.nextInt(fitToMove.size()));
			int plus = 0;
			if (Pitchforks.isPolytomy(nodeToMove) || Tools.isMultiMerger(Pitchforks.getTrueNodes(tree), nodeToMove))
				plus = 1;
			double heightMin = Math.max(nodeToMove.getChild(0).getHeight(), nodeToMove.getChild(1).getHeight());

			nodeToMove.setHeight(heightMin
					+ (chosenRecipient.getParent().getHeight() - heightMin) * Randomizer.nextDouble());


			List<Node> fitToMoveUp = getfitNodesUp(eventList, chosenRecipient.getParent().getHeight());


			logHR -= Math.log(1.0 / (chosenRecipient.getParent().getHeight() - heightMin));
			logHR -= Math.log(1.0 / fitToMove.size());
			logHR += Math.log(1.0 / (fitToMoveUp.size() + plus));

        }


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

	private List<Node> getfitNodesUp(List<GeneTreeEvent> events,
			double trNodeParentHeight) {
		
		List<Node> fitToMove = new ArrayList<>();

		for (GeneTreeEvent e : events) {
			// node is bifurcation parent, which is below transmission event and has parent
			// which is above transmission event
			if (e.type == GeneTreeEventType.BIFURCATION && e.node.getParent().getHeight() >= trNodeParentHeight
					&& e.node.getHeight() < trNodeParentHeight) {
				fitToMove.add(e.node);
			}

		}
		return fitToMove;
	}

	private List<Node> getFitNodesDown(List<GeneTreeEvent> events, double trNodeParentHeight) {
		List<Node> fitToMove = new ArrayList<>();
		for (GeneTreeEvent e : events) {
			if (e.type != GeneTreeEventType.MOCK && e.node.getHeight() == trNodeParentHeight) {
				for (int nr : e.nodesInEventNr) {
					Node parent = tree.getNode(nr);
					if (Pitchforks.isLogicalNode(parent.getChild(0)) && Pitchforks.isLogicalNode(parent.getChild(1))) {
						fitToMove.add(parent);
					}
				}
			}
		}
		return fitToMove;
	}


}
