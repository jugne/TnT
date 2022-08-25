package tnt.operators;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.SetMultimap;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import pitchfork.Pitchforks;
import starbeast2.CoordinatedOperator;
import starbeast2.SpeciesTreeInterface;
import tnt.distribution.GeneTreeEvent;
import tnt.distribution.GeneTreeIntervals;
import tnt.util.MinimumDouble;
import tnt.util.Tools;

/**
 * @author originally by Huw Ogilvie, adapted for TnT by Ugne Stolz
 */

@Description("Adjusted to also reheight the geneTree nodes which are at transmission events properly. "
		+ "Original starBeast2 description:Implements a version of the co-ordinated species and gene tree operator described in Jones (2015)."
        + "Specifically, this operator moves the species tree root node and a set of gene tree nodes related to the"
        + "species tree node by a uniform amount chosen from an exponential distribution, offset to preserve the"
        + "topology of all trees. See http://dx.doi.org/10.1101/010199 for full details.")
public class CoordinatedExponential extends CoordinatedOperator {
    public final Input<Double> betaInput = new Input<>("beta", "Beta parameter of the exponential proposal distribution", 1.0);
    public final Input<Boolean> optimiseInput = new Input<>("optimise", "Adjust beta parameter during the MCMC run to improve mixing.", true);
	public Input<List<GeneTreeIntervals>> geneTreeIntervalsInput = new Input<>("geneTreeIntervals",
			"intervals for a gene tree", new ArrayList<>(), Validate.REQUIRED);

    // scaled so that the median of the proposal distribution is equal to the mean waiting time
    // so half the time proposals will be above expectation, and half below
    private final double waitingTimeScale = 1.4426950408889634;
    protected boolean optimise;
    private double beta;
    private double lambda;
    private double waitingTime;
    private enum descendsThrough {
       LEFT_ONLY, RIGHT_ONLY, BOTH, NEITHER
    }

	SpeciesTreeInterface speciesTree;

    @Override
    public void initAndValidate() {
        beta = betaInput.get();
        lambda = 1.0 / beta;
        optimise = optimiseInput.get();
        
        speciesTree = speciesTreeInput.get();
        super.initAndValidate();
    }

    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
        // always operate on the root node
        final Node speciesTreeRoot = speciesTree.getRoot();
		// don't bother if root node is a sampled ancestor
		// TODO make it work
		if (speciesTreeRoot.isFake())
			return Double.NEGATIVE_INFINITY;

		final double currentRootHeight = speciesTreeRoot.getHeight();
		final double leftChildHeight = speciesTreeRoot.getLeft().getHeight();
		final double rightChildHeight = speciesTreeRoot.getRight().getHeight();

		final int recipientChildNr = speciesTreeRoot.getChild(1).getNr();

		final List<GeneTreeIntervals> intervalsLis = geneTreeIntervalsInput.get();
		for (GeneTreeIntervals intervals : intervalsLis) {
			HashMap<Integer, List<GeneTreeEvent>> eventList = intervals.getGeneTreeEventList();
			List<GeneTreeEvent> eventsPerTrNode = eventList.get(recipientChildNr);
			GeneTreeEvent lastEvent = eventsPerTrNode.get(eventsPerTrNode.size() - 1);

			if (currentRootHeight == lastEvent.time) {
				return Double.NEGATIVE_INFINITY;
			}
		}
        
        final MinimumDouble tipwardFreedom = new MinimumDouble();
        final SetMultimap<Integer, Node> connectingNodes = getConnectingNodes(speciesTreeRoot, tipwardFreedom);

		if (connectingNodes.isEmpty())
			return Double.NEGATIVE_INFINITY;
        tipwardFreedom.set(currentRootHeight - leftChildHeight);
        tipwardFreedom.set(currentRootHeight - rightChildHeight);
        waitingTime = tipwardFreedom.get() * waitingTimeScale;

        // the youngest age the species tree root node can be (preserving topologies)
        final double uniformShift = Randomizer.nextExponential(lambda) - tipwardFreedom.get();

        speciesTreeRoot.setHeight(currentRootHeight + uniformShift);
		Set<Node> group = new HashSet<>();
		for (Node geneTreeNode : connectingNodes.values()) {
			group.add(Pitchforks.getLogicalNode(geneTreeNode));
			group.addAll(Pitchforks.getGroup(Pitchforks.getLogicalNode(geneTreeNode)));
		}

//		int p_bf = getPolytomyCount(geneTreeInput.get().get(0));
//		int m_bf = updateMultiMergeCount(geneTreeInput.get().get(0));
//		final Tree treebf = geneTreeInput.get().get(0).copy();


		for (Node n : group) {
			if (!connectingNodes.values().contains(n)) {
				return Double.NEGATIVE_INFINITY;
			}
			n.setHeight(n.getHeight() + uniformShift);
		}

		for (Node geneTreeNode : connectingNodes.values()) {
			if (!group.contains(geneTreeNode)) {
				geneTreeNode.setHeight(geneTreeNode.getHeight() + uniformShift);
			}
		}

//		int p_af = getPolytomyCount(geneTreeInput.get().get(0));
//		int m_af = updateMultiMergeCount(geneTreeInput.get().get(0));
//		final Tree treeaf = geneTreeInput.get().get(0);
//
//		if (p_bf != p_af) {
//			System.out.println();
//		}
//		if (m_bf != m_af) {
//			return Double.NEGATIVE_INFINITY;
//		}

        // the log ratio of the density of the proposed over the current species tree root heights
        final double fLogHastingsRatio = lambda * uniformShift;

        return fLogHastingsRatio;
    }

//	private int updateMultiMergeCount(Tree tree) {
//		// Zero entries
//		int nMultiMerge = 0;
//
//		// Compute histogram
//		List<Node> trueNodes = Pitchforks.getTrueInternalNodes(tree);
//		HashMap<Double, Integer> mergerMap = new HashMap<Double, Integer>();
//
//		for (Node node : trueNodes) {
//			if (!mergerMap.keySet().contains(node.getHeight())) {
//				mergerMap.put(node.getHeight(), 0);
//			} else
//				mergerMap.put(node.getHeight(), mergerMap.get(node.getHeight()) + 1);
//		}
//
//		for (Double key : mergerMap.keySet()) {
//			if (mergerMap.get(key) > 0)
//				nMultiMerge += 1;
//		}
//		return nMultiMerge;
//	}
//
//	private int getPolytomyCount(Tree tree) {
//		int count = 0;
//
//		List<Node> trueNodes = new ArrayList<>();
//		for (Node node : tree.getNodesAsArray())
//			if (node.isRoot() || node.getParent().getHeight() > node.getHeight())
//				trueNodes.add(node);
//
//		for (Node node : trueNodes) {
//			if (!node.isLeaf() && (node.getChildren().get(0).getHeight() == node.getHeight()
//					|| node.getChildren().get(1).getHeight() == node.getHeight()))
//				count += 1;
//		}
//
//		return count;
//	}


    // identify gene tree nodes which descend through both (and also descend exclusively through)
    // the left and right children of the species tree node of interest
    private SetMultimap<Integer, Node> getConnectingNodes(Node speciesTreeNode, MinimumDouble tipwardFreedom) {
        final Node leftChildNode = speciesTreeNode.getLeft();
        final Node rightChildNode = speciesTreeNode.getRight();
        final int leftChildNodeNumber = leftChildNode.getNr();
        final int rightChildNodeNumber = rightChildNode.getNr();
        final Set<String> leftChildDescendants = findDescendants(leftChildNode, leftChildNodeNumber);
        final Set<String> rightChildDescendants = findDescendants(rightChildNode, rightChildNodeNumber);

        final SetMultimap<Integer, Node> allConnectingNodes = HashMultimap.create();
        final List<Tree> geneTrees = geneTreeInput.get();
		final List<Double> trHeights = Tools.getTransmissionHeights(speciesTree);
        for (int j = 0; j < nGeneTrees; j++) {
            final Tree geneTree = geneTrees.get(j);
			final List<Node> logical = Pitchforks.getTrueNodes(geneTree);
            final Node geneTreeRootNode = geneTree.getRoot();
            final Set<Node> jConnectingNodes = new HashSet<Node>();
            findConnectingNodes(geneTreeRootNode, jConnectingNodes, leftChildDescendants, rightChildDescendants, tipwardFreedom);
			// this for loop is the only change necessary for TNT

			for (Node n : geneTree.getNodesAsArray()) {
				if (jConnectingNodes.contains(n) && Tools.containsDoubleWithPrecision(trHeights, n.getHeight())
						&& !Tools.equalHeightWithPrecision(n, speciesTreeNode))
					return HashMultimap.create();
//					jConnectingNodes.remove(n);
				else if (!n.isLeaf() && Tools.equalHeightWithPrecision(n, speciesTreeNode)
						&& !jConnectingNodes.contains(n)) {
					return HashMultimap.create();
//					jConnectingNodes.add(n);
				} else if (jConnectingNodes.contains(n) && Tools.isMultiMerger(logical, n)) {
					return HashMultimap.create();
				}
			}

            allConnectingNodes.putAll(j, jConnectingNodes);
            geneTree.startEditing(null); // hack to stop beast.core.State.Trie memory leak
        }

        return allConnectingNodes;
    }

    private descendsThrough findConnectingNodes(Node geneTreeNode, Set<Node> connectingNodes, Set<String> leftChildDescendants, Set<String> rightChildDescendants, MinimumDouble tipwardFreedom) {
        if (geneTreeNode.isLeaf()) {
            final String descendantName = geneTreeNode.getID();
            if (leftChildDescendants.contains(descendantName)) {
                return descendsThrough.LEFT_ONLY;
            } else if (rightChildDescendants.contains(descendantName)) {
                return descendsThrough.RIGHT_ONLY;
            } else {
                return descendsThrough.NEITHER;
            }
        }

        final Node leftChild = geneTreeNode.getLeft();
        final Node rightChild = geneTreeNode.getRight();
        final descendsThrough leftDescent = findConnectingNodes(leftChild, connectingNodes, leftChildDescendants, rightChildDescendants, tipwardFreedom);
        final descendsThrough rightDescent = findConnectingNodes(rightChild, connectingNodes, leftChildDescendants, rightChildDescendants, tipwardFreedom);

        if (leftDescent == rightDescent) {
            if (leftDescent == descendsThrough.BOTH) {
                connectingNodes.add(geneTreeNode);
            }

            return leftDescent;
        }

        // this code only executes when the left and right gene tree child nodes descend through different species tree node of interest children
        final double geneTreeNodeHeight = geneTreeNode.getHeight();
        if (leftDescent == descendsThrough.BOTH) { // the gene tree node left child is a member of a connected component
            if (rightDescent == descendsThrough.NEITHER) { // the gene tree node left child is the root node of a connected component
                return descendsThrough.NEITHER;
            } else { // the gene tree node right child descends exclusively through the left XOR right child of the species tree node of interest
                // so the current gene tree node is part of a connected component but the right child is not
                final double connectedComponentTipFreedom = geneTreeNodeHeight - rightChild.getHeight();
				if (connectedComponentTipFreedom != 0)
					tipwardFreedom.set(connectedComponentTipFreedom);
                connectingNodes.add(geneTreeNode);
                return descendsThrough.BOTH;
            }
        } else if (rightDescent == descendsThrough.BOTH) { // the gene tree node right child is a member of a connected component
            if (leftDescent == descendsThrough.NEITHER) { // the gene tree node right child is the root node of a connected component
                return descendsThrough.NEITHER;
            } else { // the gene tree node left child descends exclusively through the left XOR right child of the species tree node of interest
// so the current gene tree node is part of a connected component but the left child is not
                final double connectedComponentTipFreedom = geneTreeNodeHeight - leftChild.getHeight();
				if (connectedComponentTipFreedom != 0)
					tipwardFreedom.set(connectedComponentTipFreedom);
                connectingNodes.add(geneTreeNode);
                return descendsThrough.BOTH;
            }
        } else if (leftDescent == descendsThrough.NEITHER || rightDescent == descendsThrough.NEITHER) {
            return descendsThrough.NEITHER; // the current gene tree node does not descend exclusively through the species tree node of interest
        } else { // this is a tip node of a connected component
            final double leftChildBranchLength = geneTreeNodeHeight - leftChild.getHeight();
            final double rightChildBranchLength = geneTreeNodeHeight - rightChild.getHeight();
			if (leftChildBranchLength != 0)
            tipwardFreedom.set(leftChildBranchLength);
			if (rightChildBranchLength != 0)
            tipwardFreedom.set(rightChildBranchLength);
            connectingNodes.add(geneTreeNode);
            return descendsThrough.BOTH;
        }
    }

    @Override
    public double getCoercableParameterValue() {
        return beta;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        beta = value;
        lambda = 1.0 / beta;
    }

    // optimizes beta so that it converges on the mean waiting time
    // between the first (root) and second speciation events
    @Override
    public void optimize(final double logAlpha) {
        if (optimise) {
            final double count = (m_nNrRejectedForCorrection + m_nNrAcceptedForCorrection + 1.0);
            final double delta = (waitingTime - beta) / count;
			final double beta_copy = beta;
			if (beta == -delta)
				return;
            setCoercableParameterValue(beta + delta);
//			if (beta == 0)
//				System.out.println("beta=0 in CoordinatedExponential OP");
        }
    }
}
