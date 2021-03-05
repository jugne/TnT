package tnt.operators;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import tnt.distribution.GeneTreeEvent;
import tnt.distribution.GeneTreeIntervals;

/**
 * @author Alexandra Gavryushkina
 */

@Description("Implements a narrow move between trees of different dimensions (number of nodes in trees)." +
        "It takes a random sampled node which is either a leaf with the younger sibling" +
        "or a sampled internal node. In the first case, the leaf becomes a sampled internal node by replacing its " +
        "parent (the height of the leaf remains unchanged). In the second case, the sampled internal node becomes " +
        "a leaf by inserting a new parent node at a height which is uniformly chosen on the interval " +
        " between the sampled node height and its old parent height.")
public class LeafToSampledAncestorJump extends TreeOperator {

    public Input<IntegerParameter> categoriesInput = new Input<IntegerParameter>("rateCategories", "rate category per branch");

    public Input<RealParameter> rInput =
            new Input<RealParameter>("removalProbability", "The probability of an individual to be removed from the process immediately after the sampling");
	public Input<List<GeneTreeIntervals>> geneTreeIntervalsInput = new Input<>("geneTreeIntervals",
			"intervals for a gene tree", new ArrayList<>(), Validate.REQUIRED);

    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {

        double newHeight, newRange, oldRange;
        int categoryCount = 1;
        if (categoriesInput.get() != null) {

            categoryCount = categoriesInput.get().getUpper() - categoriesInput.get().getLower() +1;
        }

        Tree tree = treeInput.get();

		Integer[] fitLeafNodeNrs = getTrNodeLeafNrsNotTransmissionOnGenes(tree);
		if (fitLeafNodeNrs.length == 0)
			return Double.NEGATIVE_INFINITY;

		Node leaf = tree.getNode(fitLeafNodeNrs[Randomizer.nextInt(fitLeafNodeNrs.length)]);
        Node parent = leaf.getParent();

        if (leaf.isDirectAncestor()) {
            oldRange = 1;
            if (parent.isRoot()) {
                final double randomNumber = Randomizer.nextExponential(1);
                newHeight = parent.getHeight() + randomNumber;
                newRange = Math.exp(randomNumber);
            } else {
                newRange = parent.getParent().getHeight() - parent.getHeight();
                newHeight = parent.getHeight() + Randomizer.nextDouble() * newRange;
            }

            if (categoriesInput.get() != null) {
                int index = leaf.getNr();
                int newValue = Randomizer.nextInt(categoryCount) + categoriesInput.get().getLower(); // from 0 to n-1, n must > 0,
                categoriesInput.get().setValue(index, newValue);
            }
        } else {
            newRange = 1;
            //make sure that the branch where a new sampled node to appear is not above that sampled node
            if (getOtherChild(parent, leaf).getHeight() >= leaf.getHeight())  {
                return Double.NEGATIVE_INFINITY;
            }
            if (parent.isRoot()) {
                oldRange = Math.exp(parent.getHeight() - leaf.getHeight());
            } else {
                oldRange = parent.getParent().getHeight() - leaf.getHeight();
            }
            newHeight = leaf.getHeight();
            if  (categoriesInput.get() != null) {
                int index = leaf.getNr();
                categoriesInput.get().setValue(index, -1);
            }
        }
        parent.setHeight(newHeight);

        //make sure that either there are no direct ancestors or r<1
        if ((rInput.get() != null) && (tree.getDirectAncestorNodeCount() > 0 && rInput.get().getValue() == 1))  {
            return Double.NEGATIVE_INFINITY;
        }

        return Math.log(newRange/oldRange);
    }

	private Integer[] getTrNodeLeafNrsNotTransmissionOnGenes(Tree transmissionTree) {
		Set<Integer> fitNodeNrs = new HashSet<Integer>();
		final List<GeneTreeIntervals> intervalsLis = geneTreeIntervalsInput.get();

		for (int n = 0; n < transmissionTree.getLeafNodeCount(); n++) {
			int recipientNr = transmissionTree.getNode(n).getParent().getChild(1).getNr();
			boolean addToFit = true;

			for (GeneTreeIntervals intervals : intervalsLis) {
				HashMap<Integer, List<GeneTreeEvent>> eventList = intervals.getGeneTreeEventList();
				List<GeneTreeEvent> eventsPerTrNode = eventList.get(recipientNr);
				GeneTreeEvent lastEvent = eventsPerTrNode.get(eventsPerTrNode.size() - 1);

				if (!transmissionTree.getNode(n).getParent().isFake()
						&& transmissionTree.getNode(n).getParent().getHeight() == lastEvent.time) {
					addToFit = false;
					break;
				}
			}

			if (addToFit) {
				fitNodeNrs.add(n);
//				fitNodeNrs.add(recipientNr);
//				fitNodeNrs.add(transmissionTree.getNode(n).getParent().getChild(0).getNr());
			}

//			int recipientNr = n.getParent().getChild(1).getNr();
//			List<GeneTreeEvent> eventsPerTrNode = eventList.get(recipientNr);
//			GeneTreeEvent lastEvent = eventsPerTrNode.get(eventsPerTrNode.size() - 1);
//
//			if (n.getParent().isFake() || n.getParent().getHeight() != lastEvent.time) {
//				fitNodeNrs.add(recipientNr);
//				fitNodeNrs.add(n.getParent().getChild(0).getNr());
//			}
		}
		return fitNodeNrs.toArray(new Integer[0]);
	}
}
