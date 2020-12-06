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
import tnt.util.Tools;

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

	public Input<SpeciesTreeInterface> transmissionTreeInput = new Input<>("transmissionTree",
			"Fully labeled transmission tree on which to simulate gene trees",
			Validate.REQUIRED);

	private Tree geneTree;
	protected SpeciesTreeInterface transmissionTree;
    
	private List<GeneTreeEvent> geneTreeEventList, storedGeneTreeEventList;
	HashMap<Integer, Integer> activeLineagesPerTransmissionTreeNode;

	// eventsPerTransmissionTreeNode is null if geneTree and transmissionTree are
	// incompatible
	HashMap<Integer, List<GeneTreeEvent>> eventsPerTransmissionTreeNode, storedEventsPerTransmissionTreeNode;

	// key: gene tree node nr, value: transmission tree node nr
	HashMap<Integer, List<Integer>> logicalGeneNodesPerTransmissionNode;

	// key: gene tree node nr, value: transmission tree node nr
	HashMap<Integer, Integer> geneTreeNodeAssignment, storedGeneTreeNodeAssignment;
//	List<Node> multifurcationParents;
	public HashMap<Integer, Integer> geneTreeTipAssignment;
	public boolean eventListDirty = true;

	List<Double> trHeights, storedTrHeights;


	@Override
	public void initAndValidate() {
		transmissionTree = transmissionTreeInput.get();
		geneTree = simulatedGeneTreeInput.get();
		if (geneTree == null) {
			geneTree = geneTreeInput.get();

			// generate map of species tree tip node names to node numbers
			final Map<String, Integer> tipNumberMap = transmissionTree.getTipNumberMap();
			geneTreeTipAssignment = new HashMap<>();
			for (int i = 0; i < geneTree.getLeafNodeCount(); i++) {
				final Node geneTreeLeafNode = geneTree.getNode(i);
				final String geneTreeLeafName = geneTreeLeafNode.getID();
				final int geneTreeLeafNumber = geneTreeLeafNode.getNr();

				if (tipNumberMap.containsKey(geneTreeLeafName)) // not in BEAUTi
					geneTreeTipAssignment.put(geneTreeLeafNumber,
							tipNumberMap.get(geneTreeLeafName));
			}
		}

		storedGeneTreeEventList = new ArrayList<>();
		activeLineagesPerTransmissionTreeNode = new HashMap<>();
		eventsPerTransmissionTreeNode = new HashMap<>();
		storedEventsPerTransmissionTreeNode = new HashMap<>();
		storedTrHeights = new ArrayList<>();

	}

	/**
	 * @return (1) list of geneTree events assigned to transmission tree branches in
	 *         the map: trTreeNodeNr -> geneTreeEventList OR (2) null if geneTree is
	 *         incompatible with transmissionTree
	 */
	public HashMap<Integer, List<GeneTreeEvent>> getGeneTreeEventList() {
		update();
		return eventsPerTransmissionTreeNode;
	}

	private void update() {
		if (!eventListDirty)
			return;

		// generate map of species tree tip node names to node numbers
		final Map<String, Integer> tipNumberMap = transmissionTree.getTipNumberMap();
		geneTreeTipAssignment = new HashMap<>();
		for (int i = 0; i < geneTree.getLeafNodeCount(); i++) {
			final Node geneTreeLeafNode = geneTree.getNode(i);
			final String geneTreeLeafName = geneTreeLeafNode.getID();
			final int geneTreeLeafNumber = geneTreeLeafNode.getNr();

			if (tipNumberMap.containsKey(geneTreeLeafName)) // not in BEAUTi
				geneTreeTipAssignment.put(geneTreeLeafNumber,
						tipNumberMap.get(geneTreeLeafName));
		}

		HashMap<Integer, List<Integer>> fakeBifurcations = new HashMap<>();
		HashMap<Double, List<Integer>> nodeTime = new HashMap<>();
		geneTreeEventList = new ArrayList<GeneTreeEvent>();
		activeLineagesPerTransmissionTreeNode = new HashMap<>();
		eventsPerTransmissionTreeNode = new HashMap<>();
		logicalGeneNodesPerTransmissionNode = new HashMap<>();
		trHeights = new ArrayList<Double>();

		geneTreeNodeAssignment = new HashMap<>();
		geneTreeNodeAssignment.putAll(geneTreeTipAssignment);

		List<Node> sortedNodes = Arrays.asList(geneTree.getNodesAsArray())
				.stream()
				.sorted(Comparator.comparingDouble(e -> e.getHeight()))
				.collect(Collectors.toList());
		for (Node n : sortedNodes) {

			if (!n.isRoot() && n.getParent().getHeight() - n.getHeight() == 0) {
				Node multifurcationParent = getMultifurcationParent(n, n.getParent());
				if (fakeBifurcations.containsKey(multifurcationParent.getNr()))
					fakeBifurcations.get(multifurcationParent.getNr()).add(n.getNr());
				else {
					List<Integer> list = new ArrayList<Integer>();
					list.add(n.getNr());
					fakeBifurcations.put(multifurcationParent.getNr(), list);
				}
			}
		}
			
		if (!Tools.fillAssignmentAndCheck(transmissionTree, geneTree.getRoot(), geneTreeNodeAssignment)) {
			eventsPerTransmissionTreeNode = null;
			return;
		}
		for (Node n : geneTree.getNodesAsArray()) {
			if (!n.isRoot() && n.getHeight() > n.getParent().getHeight()) {
				eventsPerTransmissionTreeNode = null;
				return;
			}
		}

		nodeTime = new HashMap<Double, List<Integer>>();
		for (Node node : geneTree.getNodesAsArray()) {
			if (node.isRoot() || (!node.isLeaf() && node.getParent().getHeight() != node.getHeight())) {
				if (nodeTime.containsKey(node.getHeight())) {
					nodeTime.get(node.getHeight()).add(node.getNr());
				}
				else {
					List<Integer> list = new ArrayList<Integer>();
					list.add(node.getNr());
					nodeTime.put(node.getHeight(), list);
				}
			}
		}
		
		for (Double time : nodeTime.keySet()) {
			Node first = geneTree.getNode(nodeTime.get(time).get(0));

			// check that no multiple mergers were created in different transmission
			// branches
			int trNr = geneTreeNodeAssignment.get(first.getNr());
			for (int nr : nodeTime.get(time)) {
				if (geneTreeNodeAssignment.get(nr) != trNr && (!first.isLeaf() && !geneTree.getNode(nr).isLeaf())) {
					eventsPerTransmissionTreeNode = null;
					return;
				}

			}

			GeneTreeEvent event = new GeneTreeEvent();
			event.nodesInEventNr = new ArrayList<>(nodeTime.get(time));
			event.time = first.getHeight();
			event.node = first;
			event.fakeBifCount = 0;
			if (fakeBifurcations.containsKey(first.getNr())) {
				event.nodesInEventNr.addAll(fakeBifurcations.get(first.getNr()));
				event.fakeBifCount += fakeBifurcations.get(first.getNr()).size();
				event.multiCoalCount += 1;
				event.multiCoalSize.add(fakeBifurcations.get(first.getNr()).size() + 2);
				event.type = GeneTreeEvent.GeneTreeEventType.MULTIFURCATION;
			} else {
				event.multiCoalSize.add(2);
				event.multiCoalCount += 1;
			}

			if (nodeTime.get(time).size() > 1) {

				for (Integer other : nodeTime.get(time)) {
					if (first.getNr() != other &&
								geneTreeNodeAssignment.get(first.getNr()) == geneTreeNodeAssignment
									.get(other)) {
						if (fakeBifurcations.containsKey(other)) {
							event.fakeBifCount += fakeBifurcations.get(other).size();
							event.multiCoalSize.add(fakeBifurcations.get(other).size() + 2);
							} else {
								event.multiCoalSize.add(2);
							}
							event.multiCoalCount += 1;
						event.type = GeneTreeEvent.GeneTreeEventType.MULTIFURCATION;
					} else if (geneTreeNodeAssignment.get(first.getNr()) != geneTreeNodeAssignment
							.get(other)) {
						eventsPerTransmissionTreeNode = null;
						return;
					}

				}

			} else if (event.type == null && first.getChildCount() == 2)
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





		geneTreeEventList = geneTreeEventList.stream()
				.sorted(Comparator.comparingDouble(e -> e.time))
				.collect(Collectors.toList());

        for (GeneTreeEvent event : geneTreeEventList) {

			Node trNode = transmissionTreeInput.get().getNode(geneTreeNodeAssignment.get(event.node.getNr()));
			int nrLineage = activeLineagesPerTransmissionTreeNode.get(trNode.getNr()) == null ? 0
					: activeLineagesPerTransmissionTreeNode.get(trNode.getNr());
			if (!trNode.isLeaf() && nrLineage == 0) {
				activeLineagesPerTransmissionTreeNode.put(trNode.getNr(), getLineagesRecurse(trNode));
				nrLineage = activeLineagesPerTransmissionTreeNode.get(trNode.getNr());
			}
			try {
        	switch(event.type) {
			case SAMPLE:
				nrLineage += 1;
				activeLineagesPerTransmissionTreeNode.put(trNode.getNr(), nrLineage);

				break;
			case BIFURCATION:
				nrLineage -= 1;
				activeLineagesPerTransmissionTreeNode.put(trNode.getNr(), nrLineage);
				break;
			case MULTIFURCATION:
				nrLineage -= (event.fakeBifCount + event.multiCoalCount);
				activeLineagesPerTransmissionTreeNode.put(trNode.getNr(), nrLineage);
				break;
        	}}
			catch(Exception e) {
				System.exit(-1);

			}


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
				event.type = GeneTreeEvent.GeneTreeEventType.MOCK;
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
				trHeights.add(trNode.getHeight());
			}

		eventListDirty = false;
	}

	int getLineagesRecurse(Node trNode) {
		int nLineages = 0;
		try {
		if (activeLineagesPerTransmissionTreeNode.get(trNode.getChild(0).getNr()) == null)
			nLineages += getLineagesRecurse(trNode.getChild(0));
		else
			nLineages += activeLineagesPerTransmissionTreeNode.get(trNode.getChild(0).getNr());
		if (trNode.getChildCount() > 1 && activeLineagesPerTransmissionTreeNode.get(trNode.getChild(1).getNr()) == null)
			nLineages += getLineagesRecurse(trNode.getChild(1));
		else if (trNode.getChildCount() > 1)
			nLineages += activeLineagesPerTransmissionTreeNode.get(trNode.getChild(1).getNr());
		} catch (Exception e) {
			System.exit(-1);
		}

		return nLineages;
	}

	Boolean fillAssignmentAndCheck(Node subRoot, HashMap<Integer, Integer> geneTreeNodeAssignment) {
		if (!subRoot.isLeaf()) {
			if (!fillAssignmentAndCheck(subRoot.getChild(0), geneTreeNodeAssignment))
				return false;
			Node tr1 = transmissionTree
					.getNode(geneTreeNodeAssignment.get(subRoot.getChild(0).getNr()));
			while (!tr1.isRoot()) {
				Node tr1ParentNode = tr1.getParent();
				if (tr1ParentNode.getHeight() >= subRoot.getHeight())
					break;
				tr1 = tr1ParentNode;
			}

			if (subRoot.getChildCount() > 1) {
				if (!fillAssignmentAndCheck(subRoot.getChild(1), geneTreeNodeAssignment))
					return false;
				Node tr2 = transmissionTree
						.getNode(geneTreeNodeAssignment.get(subRoot.getChild(1).getNr()));
				while (!tr2.isRoot()) {
					Node tr2ParentNode = tr2.getParent();
					if (tr2ParentNode.getHeight() >= subRoot.getHeight())
						break;
					tr2 = tr2ParentNode;
				}
				
				if (tr1.getNr() != tr2.getNr())
					return false;
			}
			if (tr1.getHeight() == subRoot.getHeight() && !subRoot.isLeaf())
				return false;
			geneTreeNodeAssignment.put(subRoot.getNr(), tr1.getNr());
		}
		return true;
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

		List<Double> tmp3 = trHeights;
		trHeights = storedTrHeights;
		storedTrHeights = tmp3;

		HashMap<Integer, Integer> tmp4 = geneTreeNodeAssignment;
		geneTreeNodeAssignment = storedGeneTreeNodeAssignment;
		storedGeneTreeNodeAssignment = tmp4;

		super.restore();
	}

	@Override
	protected void store() {
		storedGeneTreeNodeAssignment = new HashMap<Integer, Integer>(geneTreeNodeAssignment);

		storedGeneTreeEventList.clear();
		storedGeneTreeEventList.addAll(geneTreeEventList);

		storedEventsPerTransmissionTreeNode = new HashMap<Integer, List<GeneTreeEvent>>(eventsPerTransmissionTreeNode);
		storedTrHeights = new ArrayList<Double>(trHeights);
		super.store();
	}

	// DO NOT use the following without testing before!!!!
	public HashMap<Integer, Integer> getGeneTreeNodeAssignment() {
		update();
		if (eventsPerTransmissionTreeNode == null)
			return null;
		return geneTreeNodeAssignment;
	}


//
//	public HashMap<Integer, List<Integer>> getLogicalGeneNodesPerTransmissionNode() {
//		return logicalGeneNodesPerTransmissionNode;
//	}
   
}