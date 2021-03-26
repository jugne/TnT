package tnt.operators;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import pitchfork.Pitchforks;
import tnt.distribution.GeneTreeEvent;
import tnt.distribution.GeneTreeIntervals;

/**
 * @author Alexandra Gavryushkina
 */
@Description("Randomly selects true internal node (i.e. not the root and not a fake node) and move node height uniformly in interval " +
        "restricted by the node's parent and children.")
public class SAUniform extends TreeOperator {

	public Input<List<GeneTreeIntervals>> geneTreeIntervalsInput = new Input<>("geneTreeIntervals",
			"intervals for a gene tree", new ArrayList<>(), Validate.REQUIRED);

	public Input<List<Tree>> geneTreeInput = new Input<>("geneTree", "Gene tree within the species tree.",
			new ArrayList<>());

	List<GeneTreeIntervals> intervalsList;
	List<Tree> geneTrees;

    @Override
    public void initAndValidate() {
		intervalsList = geneTreeIntervalsInput.get();
		geneTrees = geneTreeInput.get();
    }

    /**
     * change the parameter and return the hastings ratio.
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
        final Tree tree = treeInput.get(this);
		intervalsList = geneTreeIntervalsInput.get();
        final int nNodeCount = tree.getNodeCount();
        int leafNodeCount = tree.getLeafNodeCount();

        //make sure that there is at least one non-fake and non-root internal node
        int fakeNodeCount = tree.getDirectAncestorNodeCount();
        if (fakeNodeCount == leafNodeCount-1 || (fakeNodeCount == leafNodeCount-2 && !tree.getRoot().isFake())) {
            return Double.NEGATIVE_INFINITY;
        }

        // randomly select internal node
        Node node;
        do {
            node = tree.getNode(leafNodeCount + Randomizer.nextInt(nNodeCount / 2));
        } while (node.isRoot() || node.isFake());

		Node recipientChild = node.getRight();
		List<Node> trueGeneNodesToMoveTogether = new ArrayList<Node>();
		double geneLower = Double.NEGATIVE_INFINITY;
		double geneUpper = Double.POSITIVE_INFINITY;
		int i = 0;
		for (GeneTreeIntervals intervals : intervalsList) {
			intervals.eventListDirty = true;
			Tree geneTree = geneTrees.get(i);
			i += 1;
			HashMap<Integer, List<GeneTreeEvent>> eventList = intervals.getGeneTreeEventList();
			List<GeneTreeEvent> eventsPerTrNode = eventList.get(recipientChild.getNr());
			GeneTreeEvent lastEvent = eventsPerTrNode.get(eventsPerTrNode.size() - 1);

			if (node.getHeight() == lastEvent.time) {

				if (lastEvent.multiCoalCount > 1) {
					for (Node n : Pitchforks.getTrueInternalNodes(geneTree)) {
						if (n.getHeight() == lastEvent.time) {
							trueGeneNodesToMoveTogether.add(n);
							double parentheight = n.getParent().getHeight();
							if (parentheight < geneUpper) {
								geneUpper = parentheight;
							}

							for (Node nChild : Pitchforks.getLogicalChildren(n)) {
								if (nChild.getHeight() > geneLower)
									geneLower = nChild.getHeight();
							}
						}
					}
				} else {
					Node logicalNode = Pitchforks.getLogicalNode(lastEvent.node);
					trueGeneNodesToMoveTogether.add(Pitchforks.getLogicalNode(lastEvent.node));
					double parentheight = logicalNode.getParent().getHeight();
					if (parentheight < geneUpper) {
						geneUpper = parentheight;
					}

					for (Node n : Pitchforks.getLogicalChildren(logicalNode)) {
						if (n.getHeight() > geneLower)
							geneLower = n.getHeight();
					}
				}
				
			}
			intervals.eventListDirty = true;
		}

		
		final double fUpper = Math.min(node.getParent().getHeight(), geneUpper);
		double fLower = Math.max(node.getLeft().getHeight(), node.getRight().getHeight());
		fLower = Math.max(fLower, geneLower);
        final double newValue = (Randomizer.nextDouble() * (fUpper - fLower)) + fLower;
        node.setHeight(newValue);

		for (Node n : trueGeneNodesToMoveTogether) {
			List<Node> group = Pitchforks.getGroup(n);
			n.setHeight(newValue);
			for (Node ng : group)
				ng.setHeight(newValue);
		}


        return 0.0;
    }

}
