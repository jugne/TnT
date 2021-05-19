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

	Integer[] geneTreeNodeAssignment, storedGeneTreeNodeAssignment;


	private Integer[] geneTreeTipAssignment;

	public HashMap<Integer, List<Integer>> inverseGeneTreeTipAssignment;
	public boolean eventListDirty = true;
	double[] trNodeOccupancy, storedTrNodeOccupancy;
	double[] trHeights, storedTrHeights;
	int nGeneNodes;
	int nTrNodes;

	boolean firstRun;


	@Override
	public void initAndValidate() {
		transmissionTree = transmissionTreeInput.get();
		geneTree = simulatedGeneTreeInput.get();
		if (geneTree == null) {
			geneTree = geneTreeInput.get();
			// generate map of species tree tip node names to node numbers
			final Map<String, Integer> tipNumberMap = transmissionTree.getTipNumberMap();
			geneTreeTipAssignment = new Integer[geneTree.getLeafNodeCount()];
			for (int i = 0; i < geneTree.getLeafNodeCount(); i++) {
				final Node geneTreeLeafNode = geneTree.getNode(i);
				final String geneTreeLeafName = geneTreeLeafNode.getID();
				final int geneTreeLeafNumber = geneTreeLeafNode.getNr();

				if (tipNumberMap.containsKey(geneTreeLeafName)) {
					geneTreeTipAssignment[geneTreeLeafNumber] = tipNumberMap.get(geneTreeLeafName);
				}

			}
		}
		
		

		storedGeneTreeEventList = new ArrayList<>();
		activeLineagesPerTransmissionTreeNode = new HashMap<>();
		eventsPerTransmissionTreeNode = new HashMap<>();
		storedEventsPerTransmissionTreeNode = new HashMap<>();
		storedTrHeights = new double[transmissionTree.getNodeCount()];
		storedGeneTreeNodeAssignment = new Integer[geneTree.getNodeCount()];

		nGeneNodes = geneTree.getNodeCount();
		nTrNodes = transmissionTree.getNodeCount();

		trNodeOccupancy = new double[nGeneNodes * nTrNodes];
		storedTrNodeOccupancy = new double[nGeneNodes * nTrNodes];

		firstRun = true;

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
		if (firstRun) {
			final Map<String, Integer> tipNumberMap = transmissionTree.getTipNumberMap();
			geneTreeTipAssignment = new Integer[geneTree.getLeafNodeCount()];
			geneTreeNodeAssignment = new Integer[geneTree.getNodeCount()];

			for (int i = 0; i < geneTree.getLeafNodeCount(); i++) {
				final Node geneTreeLeafNode = geneTree.getNode(i);
				final String geneTreeLeafName = geneTreeLeafNode.getID();


				if (tipNumberMap.containsKey(geneTreeLeafName)) {
					geneTreeTipAssignment[i] = tipNumberMap.get(geneTreeLeafName);
					geneTreeNodeAssignment[i] = tipNumberMap.get(geneTreeLeafName);
				}
			}
			trNodeOccupancy = new double[nGeneNodes * nTrNodes];
		}


		HashMap<Integer, List<Integer>> fakeBifurcations = new HashMap<>();
		HashMap<Double, List<Integer>> nodeTime = new HashMap<>();
		geneTreeEventList = new ArrayList<GeneTreeEvent>();
		activeLineagesPerTransmissionTreeNode = new HashMap<>();
		eventsPerTransmissionTreeNode = new HashMap<>();
		logicalGeneNodesPerTransmissionNode = new HashMap<>();
		trHeights = new double[transmissionTree.getNodeCount()];





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

		Boolean trTreeDirty = transmissionTree.somethingIsDirty();

		if (firstRun || trTreeDirty) {
			if (!Tools.fillAssignmentAndCheck(transmissionTree, geneTree.getRoot(), geneTreeNodeAssignment,
					trNodeOccupancy)) {
				eventsPerTransmissionTreeNode = null;
				return;
			}
		} else {
			if (!Tools.fillAndCheck(transmissionTree, geneTree.getRoot(), geneTreeNodeAssignment,
					trNodeOccupancy)) {
				eventsPerTransmissionTreeNode = null;
				return;
			}
		}

		firstRun = false;

		// DEBUG
		// Uncomment bellow for debugging

		double treeLength = getLength(geneTree);
		double treeLength2 = Arrays.stream(trNodeOccupancy).sum();
		if (Math.abs(treeLength - treeLength2) > 10e-12)
			throw new RuntimeException("lengths don't match!");

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
			int trNr = geneTreeNodeAssignment[first.getNr()];
			for (int nr : nodeTime.get(time)) {
				if (geneTreeNodeAssignment[nr] != trNr && (!first.isLeaf() && !geneTree.getNode(nr).isLeaf())) {
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
							geneTreeNodeAssignment[first.getNr()] == geneTreeNodeAssignment[other]) {
						if (fakeBifurcations.containsKey(other)) {
							event.fakeBifCount += fakeBifurcations.get(other).size();
							event.multiCoalSize.add(fakeBifurcations.get(other).size() + 2);
							} else {
								event.multiCoalSize.add(2);
							}
							event.multiCoalCount += 1;
						event.type = GeneTreeEvent.GeneTreeEventType.MULTIFURCATION;
					} else if (geneTreeNodeAssignment[first.getNr()] != geneTreeNodeAssignment[other]) {
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

			Node trNode = transmissionTreeInput.get().getNode(geneTreeNodeAssignment[event.node.getNr()]);
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
			catch (Throwable e) {
				System.out.println(e.getMessage());
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
				trHeights[trNode.getNr()] = trNode.getHeight();
			}

		eventListDirty = false;
	}

	int getLineagesRecurse(Node trNode) {
		int nLineages = 0;

		if (trNode.isLeaf() && activeLineagesPerTransmissionTreeNode.get(trNode.getNr()) == null)
			return nLineages;
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
		List<GeneTreeEvent> tmp = new ArrayList<GeneTreeEvent>(geneTreeEventList);
		geneTreeEventList = new ArrayList<GeneTreeEvent>(storedGeneTreeEventList);
		storedGeneTreeEventList = tmp;

		HashMap<Integer, List<GeneTreeEvent>> tmp2 = eventsPerTransmissionTreeNode;
		eventsPerTransmissionTreeNode = storedEventsPerTransmissionTreeNode;
		storedEventsPerTransmissionTreeNode = tmp2;

		double[] tmp3 = trHeights;
		trHeights = storedTrHeights;
		storedTrHeights = tmp3;

		Integer[] tmp4 = geneTreeNodeAssignment;
		geneTreeNodeAssignment = storedGeneTreeNodeAssignment;
		storedGeneTreeNodeAssignment = tmp4;



		double[] tmp6 = trNodeOccupancy;
		trNodeOccupancy = storedTrNodeOccupancy;
		storedTrNodeOccupancy = tmp6;

		super.restore();
	}

	@Override
	protected void store() {
		System.arraycopy(geneTreeNodeAssignment, 0, storedGeneTreeNodeAssignment, 0, geneTreeNodeAssignment.length);


		storedGeneTreeEventList.clear();
		storedGeneTreeEventList.addAll(geneTreeEventList);

		if (eventsPerTransmissionTreeNode != null)
			storedEventsPerTransmissionTreeNode = new HashMap<Integer, List<GeneTreeEvent>>(
					eventsPerTransmissionTreeNode);
		else
			storedEventsPerTransmissionTreeNode = null;


		System.arraycopy(trHeights, 0, storedTrHeights, 0, trHeights.length);
		System.arraycopy(trNodeOccupancy, 0, storedTrNodeOccupancy, 0, trNodeOccupancy.length);
		super.store();
	}

	public Integer[] getGeneTreeNodeAssignment() {
		update();
		if (eventsPerTransmissionTreeNode == null)
			return null;
		return geneTreeNodeAssignment;
	}


	public double[] getTrHeights() {
		update();
		return trHeights;
	}

	public double[] getTrNodeOccupancy() {
		update();
		return trNodeOccupancy;
	}

	private double getLength(Tree tree) {
		double length = 0;
		for (Node node : tree.getNodesAsArray()) {
			if (!node.isRoot()) {
				length += node.getLength();
			}
		}
		return length;
	}


//
//	public HashMap<Integer, List<Integer>> getLogicalGeneNodesPerTransmissionNode() {
//		return logicalGeneNodesPerTransmissionNode;
//	}
   
}