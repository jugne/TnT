package tnt.likelihood;



import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
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

	public Input<Tree> transmissionTreeInput = new Input<>("transmissionTreeInput",
			"Fully labeled transmission tree on which to simulate gene trees",
			Input.Validate.REQUIRED);

	private SimulatedGeneTree geneTree;
    
	private List<GeneTreeEvent> geneTreeEventList, storedGeneTreeEventList;
	HashMap<Node, Integer> activeLineagesPerTransmissionTreeNode;
	HashMap<Node, List<GeneTreeEvent>> eventsPerTransmissionTreeNode;

	public boolean eventListDirty = true;

	@Override
	public void initAndValidate() {
		geneTree = simulatedGeneTreeInput.get();
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

		HashMap<Double, List<Node>> fakeBifurcations = new HashMap<>();
		geneTreeEventList = new ArrayList<GeneTreeEvent>();
		activeLineagesPerTransmissionTreeNode = new HashMap<>();
		eventsPerTransmissionTreeNode = new HashMap<>();

		for (Node n : geneTree.getNodesAsArray()) {
			if (!n.isRoot() && n.getParent().getHeight() - n.getHeight() == 0) {
				if (fakeBifurcations.containsKey(n.getHeight()))
					fakeBifurcations.get(n.getHeight()).add(n);
				else {
					List<Node> list = new ArrayList<Node>();
					list.add(n);
					fakeBifurcations.put(n.getHeight(), list);
				}
				continue;
			}

			GeneTreeEvent event = new GeneTreeEvent();
			event.time = n.getHeight();
			event.node = n;

			if (n.getChildCount() == 0)
				event.type = GeneTreeEvent.GeneTreeEventType.SAMPLE;
			else if (n.getChildCount() == 2)
				event.type = GeneTreeEvent.GeneTreeEventType.BIFURCATION;
			else
				throw new RuntimeException("Network node has illegal number of children.");

			geneTreeEventList.add(event);
		}
		
		HashMap<Double, List<GeneTreeEvent>> multiCoal = new HashMap<>();

		for (GeneTreeEvent event : geneTreeEventList) {
			if (fakeBifurcations.containsKey(event.time)) {
				event.fakeBifCount = fakeBifurcations.get(event.time).size();
				event.type = GeneTreeEvent.GeneTreeEventType.MULTIFURCATION;
				event.multiCoalCount += 1;

				if (multiCoal.containsKey(event.time))
					multiCoal.get(event.time).add(event);
				else {
					List<GeneTreeEvent> coalList = new ArrayList<>();
					coalList.add(event);
					multiCoal.put(event.time, coalList);
				}
			}
		}
		
		for (double time : multiCoal.keySet()) {
			multiCoal.get(time).get(0).multiCoalSize.add(coalSize(multiCoal.get(time).get(0).node));
			if (multiCoal.get(time).size() > 1) {
				for (int i = 1; i < multiCoal.get(time).size(); i++) {
					multiCoal.get(time).get(0).multiCoalCount += multiCoal.get(time).get(i).multiCoalCount;
					multiCoal.get(time).get(0).multiCoalSize.add(coalSize(multiCoal.get(time).get(i).node));
					geneTreeEventList.remove(multiCoal.get(time).get(i));
				}
			}
		}

		int id = 0;
		for (Node trNode : transmissionTreeInput.get().getNodesAsArray()) {
			if (!trNode.isLeaf()) {
				Node mockNode = new Node();
				mockNode.setID("mock_" + id);
				id++;
				mockNode.setHeight(trNode.getHeight());
				GeneTreeEvent startEvent = new GeneTreeEvent();
				startEvent.node = mockNode;
				startEvent.type = GeneTreeEvent.GeneTreeEventType.TANSMISSION;
				startEvent.time = trNode.getHeight();
				geneTreeEventList.add(startEvent);
				geneTree.geneTreeEventAssignment.put(startEvent.node.getID(), trNode);
			}
		}

		geneTreeEventList = geneTreeEventList.stream()
				.sorted(Comparator.comparingDouble(e -> e.time))
				.collect(Collectors.toList());

        
        for (GeneTreeEvent event : geneTreeEventList) {

			Node trNode = geneTree.geneTreeEventAssignment.get(event.node.getID());
//			GeneTreeEvent startEvent = null;
			int nrLineage = activeLineagesPerTransmissionTreeNode.get(trNode) == null ? 0
					: activeLineagesPerTransmissionTreeNode.get(trNode);
			if (!trNode.isLeaf() && nrLineage == 0) {
//				if (activeLineagesPerTransmissionTreeNode.get(trNode.getChild(1)) == null)
//					System.out.println();
				activeLineagesPerTransmissionTreeNode.put(trNode, getLineagesRecurse(trNode));
				nrLineage = activeLineagesPerTransmissionTreeNode.get(trNode);
//				startEvent = new GeneTreeEvent();
//				startEvent.type = GeneTreeEvent.GeneTreeEventType.TANSMISSION;
//				startEvent.time = trNode.getHeight();
//				startEvent.lineages = nrLineage;
			}

        	switch(event.type) {
			case SAMPLE:
				nrLineage += 1;
				activeLineagesPerTransmissionTreeNode.put(trNode, nrLineage);

				break;
			case BIFURCATION:
				nrLineage -= 1;
				activeLineagesPerTransmissionTreeNode.put(trNode, nrLineage);
				if (!trNode.isRoot() && event.node.getHeight() == trNode.getParent().getHeight()) {
					Node otherChild = trNode == trNode.getParent().getChild(0) ? trNode.getParent().getChild(1)
							: trNode.getParent().getChild(0);
					activeLineagesPerTransmissionTreeNode.put(trNode.getParent(),
							nrLineage + activeLineagesPerTransmissionTreeNode.get(otherChild));
				}
				break;
			case MULTIFURCATION:
				nrLineage -= (event.fakeBifCount + event.multiCoalCount);
				activeLineagesPerTransmissionTreeNode.put(trNode, nrLineage);
				if (!trNode.isRoot() && event.node.getHeight() == trNode.getParent().getHeight()) {
					Node otherChild = trNode == trNode.getParent().getChild(0) ? trNode.getParent().getChild(1)
							: trNode.getParent().getChild(0);
					activeLineagesPerTransmissionTreeNode.put(trNode.getParent(),
							nrLineage + activeLineagesPerTransmissionTreeNode.get(otherChild));
				}
				break;
        	}
//			System.out.println("Node: " + event.node.getID());
//			System.out.println("Type: " + event.type);
//			System.out.println("Lineages after: " + lineages);
//
			event.lineages = nrLineage;

			if (eventsPerTransmissionTreeNode.get(trNode) == null) {
				List<GeneTreeEvent> l = new ArrayList<>();
				l.add(event);
				eventsPerTransmissionTreeNode.put(trNode, l);
			} else {
				eventsPerTransmissionTreeNode.get(trNode).add(event);
			}
//			if (startEvent != null) {
//				eventsPerTransmissionTreeNode.get(trNode).add(startEvent);
//			}
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

	protected int coalSize(Node n) {
		int size = 0;

		for (Node child : n.getChildren()) {
			size += countZeroBranches(child, n.getHeight());
//			if (!child.isLeaf() && !child.isFake() && child.getHeight() == n.getHeight()) {
//				size += 1;
//			}
		}

		// if multifurcation has k fake coalescent events
		// there are k+2 lineages multifurcating
		return size + 2;
	}

	protected int countZeroBranches(Node n, double height) {
		int size = 0;
		if (!n.isLeaf() && !n.isFake() && n.getHeight() == height) {
			size += countZeroBranches(n.getChild(0), height);
			size += countZeroBranches(n.getChild(1), height);
			size += 1;
		}

		return size;
	}
   
}