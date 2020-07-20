package tnt.distribution;



import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import starbeast2.SpeciesTreeInterface;
import tnt.simulator.SimulatedGeneTree;

/**
 *
 * @author Ugne Stolz
 * 
 */
@Description("Extracts the intervals from a gene tree. Points in the intervals " +
        "are defined by the heights of nodes in the tree.")
public class GeneTreeIntervals extends CalculationNode {
	public Input<Tree> geneTreeInput = new Input<Tree>("geneTree", "Gene tree for which to calculate the intervals");

	public Input<SimulatedGeneTree> simulatedGeneTreeInput=new Input<SimulatedGeneTree>("simulatedGeneTree","Gene tree for which to calculate the intervals",
			Validate.XOR, geneTreeInput);

	public Input<SpeciesTreeInterface> transmissionTreeInput = new Input<>("transmissionTreeInput",
			"Fully labeled transmission tree on which to simulate gene trees",
			Validate.REQUIRED);

	private Tree geneTree;
	protected SpeciesTreeInterface transmissionTree;
    
	private List<GeneTreeEvent> geneTreeEventList, storedGeneTreeEventList;
	HashMap<Node, Integer> activeLineagesPerTransmissionTreeNode;
	HashMap<Node, List<GeneTreeEvent>> eventsPerTransmissionTreeNode;

	HashMap<Node, List<String>> nodesPerTransmissionTreeNode;

	HashMap<String, Node> geneTreeNodeAssignment;
	public boolean eventListDirty = true;

//	HashMap<String, Node> localTipNumberMap;
//	private int[] localTipNumberMap;

	@Override
	public void initAndValidate() {
		transmissionTree = transmissionTreeInput.get();
		geneTree = simulatedGeneTreeInput.get();
		if (geneTree == null) {
			geneTree = geneTreeInput.get();

			// generate map of species tree tip node names to node numbers
			final Map<String, Integer> tipNumberMap = transmissionTree.getTipNumberMap();
			geneTreeNodeAssignment = new HashMap<>();
//			localTipNumberMap = new int[geneTree.getLeafNodeCount()];
			for (int i = 0; i < geneTree.getLeafNodeCount(); i++) {
				final Node geneTreeLeafNode = geneTree.getNode(i);
				final String geneTreeLeafName = geneTreeLeafNode.getID();
//				final int geneTreeLeafNumber = geneTreeLeafNode.getNr();

				if (tipNumberMap.containsKey(geneTreeLeafName)) // not in BEAUTi
					geneTreeNodeAssignment.put(geneTreeLeafName, transmissionTree.getNode(tipNumberMap.get(geneTreeLeafName)));
//					localTipNumberMap[geneTreeLeafNumber] = tipNumberMap.get(geneTreeLeafName);
			}
		} else {
			geneTreeNodeAssignment = simulatedGeneTreeInput.get().geneTreeSampleAssignment;
		}

		storedGeneTreeEventList = new ArrayList<>();
		activeLineagesPerTransmissionTreeNode = new HashMap<>();
		eventsPerTransmissionTreeNode = new HashMap<>();

	}

	HashMap<Node, List<GeneTreeEvent>> getGeneTreeEventList() {
		update();
		return eventsPerTransmissionTreeNode;
	}

	private void update() {
        if (!eventListDirty)
            return;

//		HashMap<Double, List<Node>> fakeBifurcations = new HashMap<>();
		HashMap<Node, List<Node>> fakeBifurcations = new HashMap<>();
		HashMap<Double, List<Node>> nodeTime = new HashMap<>();
		geneTreeEventList = new ArrayList<GeneTreeEvent>();
		activeLineagesPerTransmissionTreeNode = new HashMap<>();
		eventsPerTransmissionTreeNode = new HashMap<>();

//		int id = 0;
//		for (Node trNode : transmissionTreeInput.get().getNodesAsArray()) {
//			if (!trNode.isLeaf()) {
//				Node mockNode = new Node();
//				mockNode.setID("mock_" + id);
//				id++;
//				mockNode.setHeight(trNode.getHeight());
//				GeneTreeEvent startEvent = new GeneTreeEvent();
//				startEvent.node = mockNode;
//				startEvent.type = GeneTreeEvent.GeneTreeEventType.TANSMISSION;
//				startEvent.time = trNode.getHeight();
//				geneTreeEventList.add(startEvent);
//				geneTreeNodeAssignment.put(startEvent.node.getID(), trNode);
////				geneTree.geneTreeEventAssignment.put(startEvent.node.getID(), trNode);
//			}
//		}

		for (Node n : geneTree.getNodesAsArray()) {

			if (!n.isRoot() && n.getParent().getHeight() - n.getHeight() == 0) {
				Node multifurcationParent = getMultifurcationParent(n, n.getParent());
				if (fakeBifurcations.containsKey(multifurcationParent))
					fakeBifurcations.get(multifurcationParent).add(n);
				else {
					List<Node> list = new ArrayList<Node>();
					list.add(n);
					fakeBifurcations.put(multifurcationParent, list);
				}
			}

			Node trNode = null;
			if (!n.isLeaf()) {
				Node child = n.getChild(0);
				trNode = geneTreeNodeAssignment.get(child.getID());
				while (trNode == null) {
					child = child.getChild(0);
					trNode = geneTreeNodeAssignment.get(child.getID());
				}
				while (!trNode.isRoot() && n.getHeight() > trNode.getParent().getHeight()) {
					trNode = trNode.getParent();
				}

				// check if transmission tree is compatible with the gene tree
				// TODO make sure this works as intended
				if (n.getChildCount() > 1) {
					Node tr1 = geneTreeNodeAssignment.get(n.getChild(0).getID());
					Node tr2 = geneTreeNodeAssignment.get(n.getChild(1).getID());
					if (tr1 != tr2 && (trNode == tr1 || trNode == tr2 || !trNode.getAllChildNodesAndSelf().contains(tr1)
							|| !trNode.getAllChildNodesAndSelf().contains(tr2))) {
						geneTreeNodeAssignment = null;
						return;
					}
				}
				geneTreeNodeAssignment.put(n.getID(), trNode);
			}
		}

		nodeTime = new HashMap<Double, List<Node>>();
		for (Node node : geneTree.getNodesAsArray()) {
			if (node.isRoot() || (!node.isLeaf() && node.getParent().getHeight() != node.getHeight())) {
				if (nodeTime.containsKey(node.getHeight())) {
					nodeTime.get(node.getHeight()).add(node);
				}
				else {
					List<Node> list = new ArrayList<Node>();
					list.add(node);
					nodeTime.put(node.getHeight(), list);
				}
			}
		}
		
		for (Double time : nodeTime.keySet()) {
			Node first = nodeTime.get(time).get(0);

			GeneTreeEvent event = new GeneTreeEvent();
			event.time = first.getHeight();
			event.node = first;
			event.fakeBifCount = 0;
			if (fakeBifurcations.containsKey(first)) {
				event.fakeBifCount += fakeBifurcations.get(first).size();
				event.multiCoalCount += 1;
				event.multiCoalSize.add(fakeBifurcations.get(first).size() + 2);
				event.type = GeneTreeEvent.GeneTreeEventType.MULTIFURCATION;
//			}

			if (nodeTime.get(time).size() > 1) {

				for (Node rest : nodeTime.get(time)) {
					if (first != rest &&
							geneTreeNodeAssignment.get(first.getID()) == geneTreeNodeAssignment.get(rest.getID())) {
							if (fakeBifurcations.containsKey(rest)) {
								event.fakeBifCount += fakeBifurcations.get(rest).size();
								event.multiCoalSize.add(fakeBifurcations.get(rest).size() + 2);
							} else {
								event.multiCoalSize.add(2);
							}
							event.multiCoalCount += 1;
						event.type = GeneTreeEvent.GeneTreeEventType.MULTIFURCATION;
					}
				}
				}
			} else if (first.getChildCount() == 2)
				event.type = GeneTreeEvent.GeneTreeEventType.BIFURCATION;

			geneTreeEventList.add(event);
		}

		for (Node n : geneTree.getNodesAsArray()) {
			if (n.getChildCount() == 0) {
				GeneTreeEvent event = new GeneTreeEvent();
				event.time = n.getHeight();
				event.node = n;
				event.type = GeneTreeEvent.GeneTreeEventType.SAMPLE;
				geneTreeEventList.add(event);
			}
		}




//		for (Node n : geneTree.getNodesAsArray()) {
//			if (!n.isRoot() && n.getParent().getHeight() - n.getHeight() == 0) {
//				if (fakeBifurcations.containsKey(n.getHeight()))
//					fakeBifurcations.get(n.getHeight()).add(n);
//				else {
//					List<Node> list = new ArrayList<Node>();
//					list.add(n);
//					fakeBifurcations.put(n.getHeight(), list);
//				}
//				continue;
//			}
//
//			GeneTreeEvent event = new GeneTreeEvent();
//			event.time = n.getHeight();
//			event.node = n;
//
//			if (n.getChildCount() == 0)
//				event.type = GeneTreeEvent.GeneTreeEventType.SAMPLE;
//			else if (n.getChildCount() == 2)
//				event.type = GeneTreeEvent.GeneTreeEventType.BIFURCATION;
//			else
//				throw new RuntimeException("Network node has illegal number of children.");
//
//			geneTreeEventList.add(event);
//		}
//		
//		HashMap<Double, List<GeneTreeEvent>> multiCoal = new HashMap<>();
//
//		for (GeneTreeEvent event : geneTreeEventList) {
//			if (fakeBifurcations.containsKey(event.time)) {
//				event.fakeBifCount = fakeBifurcations.get(event.time).size();
//				event.type = GeneTreeEvent.GeneTreeEventType.MULTIFURCATION;
//				event.multiCoalCount += 1;
//
//				if (multiCoal.containsKey(event.time))
//					multiCoal.get(event.time).add(event);
//				else {
//					List<GeneTreeEvent> coalList = new ArrayList<>();
//					coalList.add(event);
//					multiCoal.put(event.time, coalList);
//				}
//			}
//		}
//		
//		for (double time : multiCoal.keySet()) {
//			multiCoal.get(time).get(0).multiCoalSize.add(coalSize(multiCoal.get(time).get(0).node));
//			if (multiCoal.get(time).size() > 1) {
//				for (int i = 1; i < multiCoal.get(time).size(); i++) {
//					multiCoal.get(time).get(0).multiCoalCount += multiCoal.get(time).get(i).multiCoalCount;
//					multiCoal.get(time).get(0).multiCoalSize.add(coalSize(multiCoal.get(time).get(i).node));
//					geneTreeEventList.remove(multiCoal.get(time).get(i));
//				}
//			}
//		}



		geneTreeEventList = geneTreeEventList.stream()
				.sorted(Comparator.comparingDouble(e -> e.time))
				.collect(Collectors.toList());

//		System.out.println(this.geneTree.getRoot().toNewick());
        for (GeneTreeEvent event : geneTreeEventList) {

//			Node trNode = null;
//        	if (!event.node.isLeaf()) {
//				Node child = event.node.getChild(0);
//				trNode = geneTreeNodeAssignment.get(child.getID());
//				while (trNode == null) {
//					child = child.getChild(0);
//					trNode = geneTreeNodeAssignment.get(child.getID());
//				}
//				while (!trNode.isRoot() && event.node.getHeight() > trNode.getParent().getHeight()) {
//					trNode = trNode.getParent();
//				}
//
//				geneTreeNodeAssignment.put(event.node.getID(), trNode);
//			} else {
			Node trNode = geneTreeNodeAssignment.get(event.node.getID());
//			}

//			Node trNode = geneTree.geneTreeEventAssignment.get(event.node.getID());
			if (trNode.getID() == null || !trNode.getID().equals("t7_1"))
				continue;
			int nrLineage = activeLineagesPerTransmissionTreeNode.get(trNode) == null ? 0
					: activeLineagesPerTransmissionTreeNode.get(trNode);
			if (!trNode.isLeaf() && nrLineage == 0) {
				activeLineagesPerTransmissionTreeNode.put(trNode, getLineagesRecurse(trNode));
				nrLineage = activeLineagesPerTransmissionTreeNode.get(trNode);
			}

        	switch(event.type) {
			case SAMPLE:
				nrLineage += 1;
				activeLineagesPerTransmissionTreeNode.put(trNode, nrLineage);
				if (trNode.getParent().isFake())
					activeLineagesPerTransmissionTreeNode.put(trNode.getParent(), nrLineage);

				break;
			case BIFURCATION:
				nrLineage -= 1;
				activeLineagesPerTransmissionTreeNode.put(trNode, nrLineage);
//				if (!trNode.isRoot() && event.node.getHeight() == trNode.getParent().getHeight()) {
//					Node otherChild = trNode == trNode.getParent().getChild(0) ? trNode.getParent().getChild(1)
//							: trNode.getParent().getChild(0);
//					activeLineagesPerTransmissionTreeNode.put(trNode.getParent(),
//							nrLineage + activeLineagesPerTransmissionTreeNode.get(otherChild));
//				}
				break;
			case MULTIFURCATION:
				nrLineage -= (event.fakeBifCount + event.multiCoalCount);
				activeLineagesPerTransmissionTreeNode.put(trNode, nrLineage);
//				if (!trNode.isRoot() && event.node.getHeight() == trNode.getParent().getHeight()) {
//					Node otherChild = trNode == trNode.getParent().getChild(0) ? trNode.getParent().getChild(1)
//							: trNode.getParent().getChild(0);
//					activeLineagesPerTransmissionTreeNode.put(trNode.getParent(),
//							nrLineage + activeLineagesPerTransmissionTreeNode.get(otherChild));
//				}
				break;
        	}

			if (nrLineage <= 0)
				System.out.println();
			event.lineages = nrLineage;

			if (eventsPerTransmissionTreeNode.get(trNode) == null) {
				List<GeneTreeEvent> l = new ArrayList<>();
				l.add(event);
				eventsPerTransmissionTreeNode.put(trNode, l);
			} else {
				eventsPerTransmissionTreeNode.get(trNode).add(event);
			}
        }

		eventListDirty = false;
	}

	int getLineagesRecurse(Node trNode) {
		int nLineages = 0;
		if (activeLineagesPerTransmissionTreeNode.get(trNode.getChild(0)) == null)
			nLineages += getLineagesRecurse(trNode.getChild(0));
		else
			nLineages += activeLineagesPerTransmissionTreeNode.get(trNode.getChild(0));
		if (activeLineagesPerTransmissionTreeNode.get(trNode.getChild(1)) == null)
			nLineages += getLineagesRecurse(trNode.getChild(1));
		else
			nLineages += activeLineagesPerTransmissionTreeNode.get(trNode.getChild(1));

		return nLineages;
	}

	Node getMultifurcationParent(Node child, Node parent) {
		Node multParent;
		if (!parent.isRoot() && child.getHeight() == parent.getParent().getHeight())
			multParent = getMultifurcationParent(child, parent.getParent());
		else
			multParent = parent;

		return multParent;
	}

	@Override
	protected boolean requiresRecalculation() {
		eventListDirty = true;
		return true;
	}

	@Override
	protected void restore() {
		List<GeneTreeEvent> tmp = geneTreeEventList;
		geneTreeEventList = storedGeneTreeEventList;
		storedGeneTreeEventList = tmp;

		super.restore();
	}

	@Override
	protected void store() {
		storedGeneTreeEventList.clear();
		storedGeneTreeEventList.addAll(geneTreeEventList);

		super.store();
	}

   
}