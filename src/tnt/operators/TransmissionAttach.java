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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import tnt.distribution.GeneTreeEvent;
import tnt.distribution.GeneTreeEvent.GeneTreeEventType;
import tnt.distribution.GeneTreeIntervals;
import tnt.util.Tools;

public class TransmissionAttach extends TreeOperator {

	public Input<GeneTreeIntervals> geneTreeIntervalsInput = new Input<>("geneTreeIntervals",
			"intervals for a gene tree", Validate.REQUIRED);

    Tree tree;
	GeneTreeIntervals intervals;
	HashMap<Integer, List<GeneTreeEvent>> eventsPerTransmissionTreeNode;
	Double heightMin;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
		intervals = geneTreeIntervalsInput.get();

    }

    @Override
    public double proposal() {

        double logHR = 0.0;
        eventsPerTransmissionTreeNode = intervals.getGeneTreeEventList();
		// fill the recipient node list on
		List<Node> recipients = new ArrayList<>();
		Node trTreeRoot = intervals.transmissionTreeInput.get().getRoot();
		getRecipients(trTreeRoot, recipients);

		int nRecipients = recipients.size();
		Node chosenRecipient = recipients.get(Randomizer.nextInt(nRecipients));
		List<GeneTreeEvent> eventList = eventsPerTransmissionTreeNode.get(chosenRecipient.getNr());
		if (eventList.size() == 0)
			return Double.NEGATIVE_INFINITY;


		heightMin = Double.NEGATIVE_INFINITY;

        if (Randomizer.nextBoolean()) {
			// Up

//			System.out.println("Up");
			for (GeneTreeEvent e : eventList) {
				// if there is already a node at transmission, don't move up.
				// we avoid this in order to not create polytomies or multiple mergers which we
				// cannot remove correctly.
				if (e.type != GeneTreeEventType.MOCK
						&& e.node.getHeight() == chosenRecipient.getParent().getHeight())
					return Double.NEGATIVE_INFINITY;
			}


			GeneTreeEvent chosenEvent = null;
			List<GeneTreeEvent> fitToMove = new ArrayList<>();
			try {
				chosenEvent = getEventUp(eventList, fitToMove,
						chosenRecipient.getParent().getHeight());
			} catch (Exception e) {
				System.out.println(e.getMessage());
				System.exit(-1);
			}

			if (chosenEvent == null)
				return Double.NEGATIVE_INFINITY;

			// record old height
			double oldHeight = chosenEvent.node.getHeight();
			// set new height for all nodes in event at
			List<Node> trueNodes = getTrueNodes(tree);
			for (int n : chosenEvent.nodesInEventNr) {
				if (Tools.isMultiMerger(trueNodes, tree.getNode(n)))
					System.out.println();
				tree.getNode(n).setHeight(chosenRecipient.getParent().getHeight());
			}
			eventList.remove(chosenEvent);

			// calc new heightMin
			getEventDown(eventList, chosenRecipient);

			if (heightMin >= oldHeight) {
				System.err.println("sanity check failed in transmission slide operator!");
				System.exit(-1);
			}

			logHR += Math.log(1.0 / (chosenRecipient.getParent().getHeight() - heightMin));
			logHR -= Math.log(1.0 / fitToMove.size());
			



			
//			logHR -=Math.log(1.0/nEventChoices);

        } else {
			// Down

//			System.out.println("Down");
			GeneTreeEvent chosenEvent = getEventDown(eventList, chosenRecipient);

			if (chosenEvent == null)
				return Double.NEGATIVE_INFINITY;


			double newHeight = heightMin
					+ (chosenRecipient.getParent().getHeight() - heightMin) * Randomizer.nextDouble();



			// set new height for all nodes in event at
			for (int n : chosenEvent.nodesInEventNr) {
				tree.getNode(n).setHeight(newHeight);
			}
			
			List<GeneTreeEvent> fitToMove = new ArrayList<>();
			getEventUp(eventList, fitToMove,
						chosenRecipient.getParent().getHeight());

			logHR -= Math.log(1.0 / (chosenRecipient.getParent().getHeight() - heightMin));
			logHR += Math.log(1.0 / fitToMove.size());

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

	private GeneTreeEvent getEventUp(List<GeneTreeEvent> events, List<GeneTreeEvent> fitToMode,
			double trNodeParentHeight) {
		// Get first event below trNodeHeight which is not SAMPLE or mock TRANSMISSION
		// event
		heightMin = Double.NEGATIVE_INFINITY;
		GeneTreeEvent first = null;
		for (int i = events.size() - 1; i >= 0; i--) {
			if (events.get(i).type != GeneTreeEventType.SAMPLE
					&& events.get(i).type != GeneTreeEventType.MOCK) {
				if (trNodeParentHeight == events.get(i).node.getHeight()) // do not allow the move if there is already a
																			// an event at transmission height
					return null;

				if (first == null)
					first = events.get(i);
				fitToMode.add(events.get(i));

			}
		}

		return first;


//		List<GeneTreeEvent> eligible = new ArrayList<>();
//		for (GeneTreeEvent e : events) {
//			if (!e.node.isLeaf() && trNodeHeight != e.node.getHeight()) {
//				eligible.add(e);
//				nToChoose += 1;
//				
//			}
//		}
//		if (eligible.isEmpty())
//			return null;
//
//		return eligible.get(Randomizer.nextInt(nToChoose));
	}

	private GeneTreeEvent getEventDown(List<GeneTreeEvent> events, Node trNode) {
		// Get first event below trNodeHeight which is not SAMPLE or mock TRANSMISSION
		// event
		heightMin = Double.NEGATIVE_INFINITY;
		GeneTreeEvent e = null;
		for (int i = events.size() - 1; i >= 0; i--) {
			if (heightMin == Double.NEGATIVE_INFINITY && events.get(i).type != GeneTreeEventType.MOCK) {
				if (trNode.getParent().getHeight() != events.get(i).node.getHeight())
					heightMin = events.get(i).node.getHeight();
				else if (e == null && trNode.getParent().getHeight() == events.get(i).node.getHeight()) {
					e = events.get(i);
				}
				if (heightMin != Double.NEGATIVE_INFINITY && e != null)
					break;

			}
		}
		
		if (heightMin==Double.NEGATIVE_INFINITY)
			heightMin = trNode.getHeight();
		
		if (e == null)
			return null;
		return e;
	}
}
