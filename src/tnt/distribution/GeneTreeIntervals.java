package tnt.distribution;



import java.util.ArrayList;
import java.util.Arrays;
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
	HashMap<Integer, Integer> activeLineagesPerTransmissionTreeNode;
	HashMap<Integer, List<GeneTreeEvent>> eventsPerTransmissionTreeNode, storedEventsPerTransmissionTreeNode;

	// key: gene tree node nr, value: transmission tree node nr
	HashMap<Integer, List<Integer>> logicalGeneNodesPerTransmissionNode;

	// key: gene tree node nr, value: transmission tree node nr
	HashMap<Integer, Integer> geneTreeNodeAssignment, storeGeneTreeNodeAssignement;
	List<Node> multifurcationParents;
	HashMap<Integer, Integer> geneTreeTipAssignment;
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
			geneTreeTipAssignment = new HashMap<>();
//			localTipNumberMap = new int[geneTree.getLeafNodeCount()];
			for (int i = 0; i < geneTree.getLeafNodeCount(); i++) {
				final Node geneTreeLeafNode = geneTree.getNode(i);
				final String geneTreeLeafName = geneTreeLeafNode.getID();
				final int geneTreeLeafNumber = geneTreeLeafNode.getNr();

				if (tipNumberMap.containsKey(geneTreeLeafName)) // not in BEAUTi
					geneTreeTipAssignment.put(geneTreeLeafNumber,
							tipNumberMap.get(geneTreeLeafName));
//					localTipNumberMap[geneTreeLeafNumber] = tipNumberMap.get(geneTreeLeafName);
			}
//		} else {
//			geneTreeNodeAssignment = simulatedGeneTreeInput.get().geneTreeSampleAssignment;
		}

		storedGeneTreeEventList = new ArrayList<>();
		activeLineagesPerTransmissionTreeNode = new HashMap<>();
		eventsPerTransmissionTreeNode = new HashMap<>();
		storedEventsPerTransmissionTreeNode = new HashMap<>();

	}

	public HashMap<Integer, List<GeneTreeEvent>> getGeneTreeEventList() {
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
		logicalGeneNodesPerTransmissionNode = new HashMap<>();

		geneTreeNodeAssignment = new HashMap<>();
		geneTreeNodeAssignment.putAll(geneTreeTipAssignment);
		multifurcationParents = new ArrayList<Node>();
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

		List<Node> sortedNodes = Arrays.asList(geneTree.getNodesAsArray())
				.stream()
				.sorted(Comparator.comparingDouble(e -> e.getHeight()))
				.collect(Collectors.toList());
		for (Node n : sortedNodes) {

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
				// fake nodes at multifurcation can be sorted wrong. Get first real child.
				while (child.getHeight() == n.getHeight())
					child = child.getChild(0);
				trNode = transmissionTreeInput.get().getNode(geneTreeNodeAssignment.get(child.getNr()));
				while (trNode == null) {
					child = child.getChild(0);
					trNode = transmissionTreeInput.get().getNode(geneTreeNodeAssignment.get(child.getNr()));
				}
				while (!trNode.isRoot() && n.getHeight() > trNode.getParent().getHeight()) {
					trNode = trNode.getParent();
				}

				// check if transmission tree is compatible with the gene tree
				// polytomies make nodes at the same height which might not be processed in
				// child-parent order.
				// Therefore, we check them separately after node assignment is done.

				// TODO make sure this works as intended
				if (n.getChildCount() > 1 &&
						n.getChild(0).getHeight() != n.getHeight() && n.getChild(1).getHeight() != n.getHeight()) {
					Node tr1 = transmissionTreeInput.get()
							.getNode(geneTreeNodeAssignment.get(n.getChild(0).getNr()));
					Node tr2 = transmissionTreeInput.get()
							.getNode(geneTreeNodeAssignment.get(n.getChild(1).getNr()));
					if (tr1 != tr2 && ((trNode == tr1 && !trNode.getAllChildNodesAndSelf().contains(tr2))
							|| (trNode == tr2 && !trNode.getAllChildNodesAndSelf().contains(tr1))
							|| (!trNode.getAllChildNodesAndSelf().contains(tr2)
									&& !trNode.getAllChildNodesAndSelf().contains(tr1)))) {
						eventsPerTransmissionTreeNode = null;
						return;
					}
				}
				if (n.getChildCount() > 1 && (n.isRoot() || n.getParent().getHeight() != n.getHeight())
						&& (n.getChild(0).getHeight() == n.getHeight() || n.getChild(1).getHeight() == n.getHeight()))
					multifurcationParents.add(n);

				geneTreeNodeAssignment.put(n.getNr(), trNode.getNr());
			}
		}

		for (Node m : multifurcationParents) {
			Node trNode = transmissionTree.getNode(geneTreeNodeAssignment.get(m.getNr()));
			Node tr1 = transmissionTreeInput.get().getNode(geneTreeNodeAssignment.get(m.getChild(0).getNr()));
			Node tr2 = transmissionTreeInput.get().getNode(geneTreeNodeAssignment.get(m.getChild(1).getNr()));
			if (tr1 != tr2 && ((trNode == tr1 && !trNode.getAllChildNodesAndSelf().contains(tr2))
					|| (trNode == tr2 && !trNode.getAllChildNodesAndSelf().contains(tr1))
					|| (!trNode.getAllChildNodesAndSelf().contains(tr2)
							&& !trNode.getAllChildNodesAndSelf().contains(tr1)))) {
				eventsPerTransmissionTreeNode = null;
				return;
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
			event.node = first;// Pitchforks.getLogicalNode(first); // TODO validate this
			event.fakeBifCount = 0;
			if (fakeBifurcations.containsKey(first)) {
				event.fakeBifCount += fakeBifurcations.get(first).size();
				event.multiCoalCount += 1;
				event.multiCoalSize.add(fakeBifurcations.get(first).size() + 2);
				event.type = GeneTreeEvent.GeneTreeEventType.MULTIFURCATION;
			}

			if (nodeTime.get(time).size() > 1) {

				for (Node rest : nodeTime.get(time)) {
					if (first != rest &&
								geneTreeNodeAssignment.get(first.getNr()) == geneTreeNodeAssignment
										.get(rest.getNr())) {
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

//				if (event.multiCoalSize.size() > 2)
//					System.out.println();
			} else if (event.type == null && first.getChildCount() == 2)
				event.type = GeneTreeEvent.GeneTreeEventType.BIFURCATION;

//			Node trNode = geneTreeNodeAssignment.get(event.node.getNr());
//			if (!trNode.isRoot()
//					&& (trNode.getParent().getChild(0) == trNode && !trNode.getParent().isFake())
//					&& event.time == trNode.getParent().getHeight())
//				event.type = GeneTreeEvent.GeneTreeEventType.TANSMISSION;
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

//		for (Node trNode : transmissionTree.getInternalNodes()) {
//			GeneTreeEvent event = new GeneTreeEvent();
//			event.time = trNode.getHeight();
//			event.node = new Node();
//			event.node.setNr(geneTree.getNodeCount() + 100);
//			geneTreeNodeAssignment.put(event.node.getNr(), trNode);
//			event.type = GeneTreeEvent.GeneTreeEventType.TRANSMISSION;
//			geneTreeEventList.add(event);
//		}




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
			Node trNode = transmissionTreeInput.get().getNode(geneTreeNodeAssignment.get(event.node.getNr()));
//			}

//			Node trNode = geneTree.geneTreeEventAssignment.get(event.node.getID());
//			if (trNode.getID() == null || !trNode.getID().equals("t7_1"))
//				continue;
			int nrLineage = activeLineagesPerTransmissionTreeNode.get(trNode.getNr()) == null ? 0
					: activeLineagesPerTransmissionTreeNode.get(trNode.getNr());
			if (!trNode.isLeaf() && nrLineage == 0) {
				activeLineagesPerTransmissionTreeNode.put(trNode.getNr(), getLineagesRecurse(trNode));
				nrLineage = activeLineagesPerTransmissionTreeNode.get(trNode.getNr());
			}

        	switch(event.type) {
			case SAMPLE:
				nrLineage += 1;
				activeLineagesPerTransmissionTreeNode.put(trNode.getNr(), nrLineage);
//				if (trNode.getParent() != null && trNode.getParent().isFake()) // null can only happen on transmission
//																				// trees with a single branch
//					activeLineagesPerTransmissionTreeNode.put(trNode.getParent().getNr(), nrLineage);

				break;
			case BIFURCATION:
				nrLineage -= 1;
				activeLineagesPerTransmissionTreeNode.put(trNode.getNr(), nrLineage);
//				if (!trNode.isRoot() && event.node.getHeight() == trNode.getParent().getHeight()) {
//					Node otherChild = trNode == trNode.getParent().getChild(0) ? trNode.getParent().getChild(1)
//							: trNode.getParent().getChild(0);
//					activeLineagesPerTransmissionTreeNode.put(trNode.getParent(),
//							nrLineage + activeLineagesPerTransmissionTreeNode.get(otherChild));
//				}
				break;
			case MULTIFURCATION:
				nrLineage -= (event.fakeBifCount + event.multiCoalCount);
				activeLineagesPerTransmissionTreeNode.put(trNode.getNr(), nrLineage);
//				if (!trNode.isRoot() && event.node.getHeight() == trNode.getParent().getHeight()) {
//					Node otherChild = trNode == trNode.getParent().getChild(0) ? trNode.getParent().getChild(1)
//							: trNode.getParent().getChild(0);
//					activeLineagesPerTransmissionTreeNode.put(trNode.getParent(),
//							nrLineage + activeLineagesPerTransmissionTreeNode.get(otherChild));
//				}
				break;
//			case TRANSMISSION:
//				if (!trNode.isRoot() && !trNode.isFake() && !trNode.isLeaf())
//					activeLineagesPerTransmissionTreeNode.put(trNode, nrLineage);
//				break;
        	}

			if (nrLineage <= 0)
				System.out.println();
			event.lineages = nrLineage;

			if (eventsPerTransmissionTreeNode.get(trNode.getNr()) == null) {
				List<GeneTreeEvent> l = new ArrayList<>();
				l.add(event);
				eventsPerTransmissionTreeNode.put(trNode.getNr(), l);
			} else {
				eventsPerTransmissionTreeNode.get(trNode.getNr()).add(event);
			}
			
			if (logicalGeneNodesPerTransmissionNode.get(trNode.getNr()) == null) {
				List<Integer> l = new ArrayList<>();
				l.add(event.node.getNr());
				logicalGeneNodesPerTransmissionNode.put(trNode.getNr(), l);
			} else {
				logicalGeneNodesPerTransmissionNode.get(trNode.getNr()).add(event.node.getNr());
			}

        }

		for (Node trNode : transmissionTree.getNodesAsArray())
			if (!trNode.isLeaf()) { // !trNode.isFake()
				GeneTreeEvent event = new GeneTreeEvent();
				event.time = trNode.getHeight();
				event.type = GeneTreeEvent.GeneTreeEventType.TRANSMISSION;
				event.lineages += getLineagesRecurse(trNode);

				if (eventsPerTransmissionTreeNode.get(trNode.getNr()) == null) {
					List<GeneTreeEvent> l = new ArrayList<>();
					l.add(event);
					eventsPerTransmissionTreeNode.put(trNode.getNr(), l);
				} else {
					eventsPerTransmissionTreeNode.get(trNode.getNr()).add(0, event); // transmission event here is
																						// marked at the start of the
																						// branch going backwards in
																						// time.
					// It is a mock event, used only to correctly store the number of incoming gene
					// lineages from child transmission tree branches to parent branch.
				}

			}

		eventListDirty = false;
	}

	int getLineagesRecurse(Node trNode) {
		int nLineages = 0;
		if (activeLineagesPerTransmissionTreeNode.get(trNode.getChild(0).getNr()) == null)
			nLineages += getLineagesRecurse(trNode.getChild(0));
		else
			nLineages += activeLineagesPerTransmissionTreeNode.get(trNode.getChild(0).getNr());
		if (trNode.getChildCount() > 1 && activeLineagesPerTransmissionTreeNode.get(trNode.getChild(1).getNr()) == null)
			nLineages += getLineagesRecurse(trNode.getChild(1));
		else if (trNode.getChildCount() > 1)
			nLineages += activeLineagesPerTransmissionTreeNode.get(trNode.getChild(1).getNr());

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

		HashMap<Integer, List<GeneTreeEvent>> tmp2 = eventsPerTransmissionTreeNode;
		eventsPerTransmissionTreeNode = storedEventsPerTransmissionTreeNode;
		storedEventsPerTransmissionTreeNode = tmp2;

		super.restore();
	}

	@Override
	protected void store() {
		storedGeneTreeEventList.clear();
		storedGeneTreeEventList.addAll(geneTreeEventList);

		storedEventsPerTransmissionTreeNode = new HashMap<Integer, List<GeneTreeEvent>>(eventsPerTransmissionTreeNode);
		super.store();
	}

	public HashMap<Integer, Integer> getGeneTreeNodeAssignment() {
		update();
		return geneTreeNodeAssignment;
	}

	public HashMap<Integer, List<Integer>> getLogicalGeneNodesPerTransmissionNode() {
		return logicalGeneNodesPerTransmissionNode;
	}
   
}