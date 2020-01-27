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

/**
 *
 * @author Ugne Stolz
 * 
 */
@Description("Extracts the intervals from a gene tree. Points in the intervals " +
        "are defined by the heights of nodes in the tree.")
public class GeneTreeIntervals extends CalculationNode {
	public Input<Tree> geneTreeInput = new Input<Tree>("geneTree", "Gene tree for which to calculate the intervals",
			Validate.REQUIRED);

    private Tree geneTree;
    
	private List<GeneTreeEvent> geneTreeEventList, storedGeneTreeEventList;

	public boolean eventListDirty = true;

	@Override
	public void initAndValidate() {
		geneTree = geneTreeInput.get();
		storedGeneTreeEventList = new ArrayList<>();
	}

	List<GeneTreeEvent> getGeneTreeEventList() {
		update();
		return geneTreeEventList;
	}

	private void update() {
        if (!eventListDirty)
            return;

		HashMap<Double, List<Node>> fakeBifurcations = new HashMap<>();
		geneTreeEventList = new ArrayList<GeneTreeEvent>();

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
			if (multiCoal.get(time).size() > 1) {
				for (int i = 1; i < multiCoal.get(time).size(); i++) {
					multiCoal.get(time).get(0).multiCoalCount += multiCoal.get(time).get(i).multiCoalCount;
					geneTreeEventList.remove(multiCoal.get(time).get(i));
				}
			}
		}
		
		
		for (GeneTreeEvent e : geneTreeEventList) {
			e.time = Math.round(e.time * 10000000000.0) / 10000000000.0;
		}


		geneTreeEventList = geneTreeEventList.stream()
				.sorted(Comparator.comparingDouble(e -> e.time))
				.collect(Collectors.toList());


        int lineages = 0;
        
        for (GeneTreeEvent event : geneTreeEventList) {
        	switch(event.type) {
			case SAMPLE:
				lineages += 1;
				break;
			case BIFURCATION:
				lineages -= 1;
				break;
			case MULTIFURCATION:
				lineages = lineages - (event.fakeBifCount + event.multiCoalCount);
				break;
        	}
//			System.out.println("Node: " + event.node.getID());
//			System.out.println("Type: " + event.type);
//			System.out.println("Lineages after: " + lineages);
//
			if (lineages < 2)
				System.out.println(geneTree.getRoot().toNewick());
			event.lineages = lineages;
        }

		eventListDirty = false;
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