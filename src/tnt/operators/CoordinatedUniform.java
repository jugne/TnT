package tnt.operators;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.SetMultimap;

import beast.core.Description;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import pitchfork.Pitchforks;
import starbeast2.CoordinatedOperator;
import starbeast2.SpeciesTreeInterface;
import tnt.util.MinimumDouble;
import tnt.util.Tools;

/**
 * @author originally by Huw Ogilvie, adapted for TnT by Ugne Stolz
 */

@Description("Adjusted to also reheight the geneTree nodes which are at transmission events properly. "
		+ "Original starBeast2 description: Implements a version of the co-ordinated species and gene tree operator described in Jones (2015)."
        + "Specifically, this operator moves a species tree node and a set of gene tree nodes related to the"
        + "species tree node by a uniform amount chosen from a range which preserves the topology of all trees."
        + "See http://dx.doi.org/10.1101/010199 for full details.")
public class CoordinatedUniform extends CoordinatedOperator {
    private enum descendsThrough {
       LEFT_ONLY, RIGHT_ONLY, BOTH, NEITHER
    }

    SpeciesTreeInterface speciesTree;
	List<Double> connectingNodesHeights;

    @Override
    public void initAndValidate() {
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
        final double fLogHastingsRatio = 0.0; // this move is uniform in both directions


        // don't operate on sampled ancestor nodes
        // TODO make it work
        final int nLeaves = speciesTree.getLeafNodeCount();
        final int nNodes = nLeaves + nLeaves - 1;
        final Node[] nodeArray = speciesTree.getNodesAsArray();
        final Node[] trueNonRootInternalNodes = new Node[nNodes];
		connectingNodesHeights = new ArrayList<Double>();

        int trueNonRootInternalNodeCount = 0;
        for (int nodeIndex = nLeaves; nodeIndex < nNodes; nodeIndex++) {
            final Node node = nodeArray[nodeIndex];
            if (!(node.isFake() || node.isRoot())) {
                trueNonRootInternalNodes[trueNonRootInternalNodeCount] = node;
                trueNonRootInternalNodeCount++;
            }
        }

        if (trueNonRootInternalNodeCount == 0) {
            return Double.NEGATIVE_INFINITY;
        }

        final Node speciesTreeNode = trueNonRootInternalNodes[Randomizer.nextInt(trueNonRootInternalNodeCount)];
        final double speciesTreeNodeHeight = speciesTreeNode.getHeight();

        final MinimumDouble tipwardFreedom = new MinimumDouble();
        final MinimumDouble rootwardFreedom = new MinimumDouble();
        final SetMultimap<Integer, Node> connectingNodes = getConnectingNodes(speciesTreeNode, tipwardFreedom, rootwardFreedom);

		if (connectingNodes.isEmpty())
			return Double.NEGATIVE_INFINITY;

        final double leftChildBranchLength = speciesTreeNodeHeight - speciesTreeNode.getLeft().getHeight();
        final double rightChildBranchLength = speciesTreeNodeHeight - speciesTreeNode.getRight().getHeight();
        final double speciesTreeNodeBranchLength = speciesTreeNode.getParent().getHeight() - speciesTreeNodeHeight;
        tipwardFreedom.set(leftChildBranchLength);
        tipwardFreedom.set(rightChildBranchLength);
        rootwardFreedom.set(speciesTreeNodeBranchLength);

        final double twf = tipwardFreedom.get();
        final double rwf = rootwardFreedom.get();
        final double uniformShift = (Randomizer.nextDouble() * (twf + rwf)) - twf;

        speciesTreeNode.setHeight(speciesTreeNode.getHeight() + uniformShift);
		speciesTreeNode.makeDirty(Tree.IS_FILTHY);
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
				geneTreeNode.makeDirty(Tree.IS_FILTHY);
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
    private SetMultimap<Integer, Node> getConnectingNodes(Node speciesTreeNode, MinimumDouble tipwardFreedom, MinimumDouble rootwardFreedom) {
        final Node leftChildNode = speciesTreeNode.getLeft();
        final Node rightChildNode = speciesTreeNode.getRight();
        final int leftChildNodeNumber = leftChildNode.getNr();
        final int rightChildNodeNumber = rightChildNode.getNr();
		// find descendants gets only leaf decendants
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
            findConnectingNodes(geneTreeRootNode, jConnectingNodes, leftChildDescendants, rightChildDescendants, tipwardFreedom, rootwardFreedom);
			// this for loop is the only change necessary for TNT
			for (Node n : geneTree.getNodesAsArray()) {
				if (jConnectingNodes.contains(n) && Tools.containsDoubleWithPrecision(trHeights, n.getHeight())
						&& !Tools.equalHeightWithPrecisionNode(n, speciesTreeNode))
					return HashMultimap.create();
//					jConnectingNodes.remove(n);
				else if (!n.isLeaf() && Tools.equalHeightWithPrecisionNode(n, speciesTreeNode)
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

    private descendsThrough findConnectingNodes(Node geneTreeNode, Set<Node> connectingNodes, Set<String> leftChildDescendants, Set<String> rightChildDescendants, MinimumDouble tipwardFreedom, MinimumDouble rootwardFreedom) {
    	
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
        final descendsThrough leftDescent = findConnectingNodes(leftChild, connectingNodes, leftChildDescendants, rightChildDescendants, tipwardFreedom, rootwardFreedom);
        final descendsThrough rightDescent = findConnectingNodes(rightChild, connectingNodes, leftChildDescendants, rightChildDescendants, tipwardFreedom, rootwardFreedom);

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
				final double connectedComponentRootFreedom = geneTreeNodeHeight - leftChild.getHeight();
				rootwardFreedom.set(connectedComponentRootFreedom);
                return descendsThrough.NEITHER;
            } else { // the gene tree node right child descends exclusively through the left XOR right child of the species tree node of interest
                // so the current gene tree node is part of a connected component but the right child is not
                final double connectedComponentDescendantBranchLength = geneTreeNodeHeight - rightChild.getHeight();
				tipwardFreedom.set(connectedComponentDescendantBranchLength);
				connectingNodes.add(geneTreeNode);
                return descendsThrough.BOTH;
			}
        } else if (rightDescent == descendsThrough.BOTH) { // the gene tree node right child is a member of a connected component
            if (leftDescent == descendsThrough.NEITHER) { // the gene tree node right child is the root node of a connected component
				final double connectedComponentRootFreedom = geneTreeNodeHeight - rightChild.getHeight();
				rootwardFreedom.set(connectedComponentRootFreedom);
                return descendsThrough.NEITHER;
            } else { // the gene tree node left child descends exclusively through the left XOR right child of the species tree node of interest
// so the current gene tree node is part of a connected component but the left child is not
                final double connectedComponentTipFreedom = geneTreeNodeHeight - leftChild.getHeight();
				tipwardFreedom.set(connectedComponentTipFreedom);
				connectingNodes.add(geneTreeNode);
                return descendsThrough.BOTH;
            }
        } else if (leftDescent == descendsThrough.NEITHER || rightDescent == descendsThrough.NEITHER) {
            return descendsThrough.NEITHER; // the current gene tree node does not descend exclusively through the species tree node of interest
        } else { // this is a tip node of a connected component
            final double leftChildBranchLength = geneTreeNodeHeight - leftChild.getHeight();
			final double rightChildBranchLength = geneTreeNodeHeight - rightChild.getHeight();
			tipwardFreedom.set(leftChildBranchLength);
			tipwardFreedom.set(rightChildBranchLength);
			connectingNodes.add(geneTreeNode);
            return descendsThrough.BOTH;
        }
    }
}
