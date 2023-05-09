package tnt.util;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import pitchfork.Pitchforks;
import starbeast2.SpeciesTreeInterface;
import tnt.distribution.GeneTreeEvent;
import tnt.distribution.GeneTreeIntervals;

public class Tools {

	public final static double globalPrecisionThreshold = 1e-10;

	public static boolean isMultiMerger(List<Node> allLogicalNodes, Node logicalNode) {
		for (Node n : allLogicalNodes) {
			if (n.getHeight() == logicalNode.getHeight() && n.getNr() != logicalNode.getNr())
				return true;
		}
		return false;
	}

	public static List<Node> getGeneNodesAtTransmissionWithPrecision(List<Node> nodeList, List<Double> transmissionHeights) {
		List<Node> tmp = new ArrayList<>();
		for (Node n : nodeList) {
			if (containsDoubleWithPrecision(transmissionHeights, n.getHeight()))
				tmp.add(n);
		}

		return tmp;
	}


//	public static List<Node> getGeneNodesNotAtTransmission(List<Node> nodeList, List<Double> transmissionHeights) {
//		List<Node> tmp = new ArrayList<>();
//		for (Node n : nodeList) {
//			if (!containsDoubleWithPrecision(transmissionHeights, n.getHeight()))
//				tmp.add(n);
//		}
//
//		return tmp;
//	}

	public static List<Node> getGeneNodesWithParentNotAtTransmission(List<Node> nodeList,
			List<Double> transmissionHeights) {
		List<Node> tmp = new ArrayList<>();
		for (Node n : nodeList) {
			if (n.isRoot() || !containsDoubleWithPrecision(transmissionHeights, n.getParent().getHeight()))
				tmp.add(n);
		}

		return tmp;
	}

	public static Integer[] getTrNodeNrsNotTransmissionOnGenes(Tree transmissionTree,
			List<GeneTreeIntervals> geneTreeIntervals, boolean bothChildren) {
		Set<Integer> fitNodeNrs = new HashSet<Integer>();
		final List<GeneTreeIntervals> intervalsLis = geneTreeIntervals;

		for (Node n : transmissionTree.getNodesAsArray()) {
			if (n.isRoot()) {
				fitNodeNrs.add(n.getNr());
				continue;
			}

			int recipientNr = n.getParent().getChild(1).getNr();
			boolean addToFit = true;

			for (GeneTreeIntervals intervals : intervalsLis) {
				HashMap<Integer, List<GeneTreeEvent>> eventList = intervals.getGeneTreeEventList();
				List<GeneTreeEvent> eventsPerTrNode = eventList.get(recipientNr);
				if (eventsPerTrNode == null || eventsPerTrNode.size()==0)
					continue;
				GeneTreeEvent lastEvent = eventsPerTrNode.get(eventsPerTrNode.size() - 1);

				if (!n.getParent().isFake()
						&& Tools.equalWithPrecision(n.getParent().getHeight(), lastEvent.time)) {
					addToFit = false;
					break;
				}
			}

			if (addToFit) {
				fitNodeNrs.add(recipientNr);
				fitNodeNrs.add(n.getParent().getChild(0).getNr());
			} else if (!bothChildren) {// see if this works
				fitNodeNrs.add(n.getParent().getChild(0).getNr());
			}

		}

		return fitNodeNrs.toArray(new Integer[0]);
	}

	/**
	 * TRANSMISSION TREE METHOD
	 * 
	 * @param trNode Transmission tree node for which to get all upstream nodes in a
	 *               tree, starting with its parent
	 * @return list of all upstream nodes starting at trNode
	 */
	public static List<Integer> getAllParentNrs(Node trNode) {
		List<Integer> tmp = new ArrayList<>();
		Node trNodeTmp = trNode;
		while (!trNodeTmp.isRoot()) {
			tmp.add(trNodeTmp.getParent().getNr());
			trNodeTmp = trNodeTmp.getParent();
		}

		return tmp;
	}

	/**
	 * TRANSMISSION TREE METHOD
	 * 
	 * @param trNode Transmission tree node for which to get all downstream nodes in
	 *               a tree, starting with its child nodes
	 * @return List of all downstream nodes starting at trNode
	 */
	public static List<Integer> getChildNrs(Node trNode) {
		List<Integer> tmp = new ArrayList<>();

		if (trNode.isLeaf())
			return tmp;
		for (Node child : trNode.getChildren()) {
			tmp.add(child.getNr());
			tmp.addAll(getChildNrs(child));
		}

		return tmp;
	}

	public static List<Node> getGeneRootAndNodesWithParentsAtTransmission(List<Node> nodeList,
			List<Double> transmissionHeights) {
		List<Node> tmp = new ArrayList<>();
		for (Node n : nodeList) {
			if (n.isRoot())
				tmp.add(n);
			else if (containsDoubleWithPrecision(transmissionHeights, n.getParent().getHeight()))
				tmp.add(n);
		}

		return tmp;
	}

	public static <T> boolean listEqualsIgnoreOrder(List<T> list1, List<T> list2) {
		return new HashSet<>(list1).equals(new HashSet<>(list2));
	}

	public static List<Double> getTransmissionHeights(SpeciesTreeInterface transmissionTree) {
		List<Double> trHeights = new ArrayList<>();
		for (Node trNode : transmissionTree.getInternalNodes()) {
			if (!trNode.isFake()) {
				trHeights.add(trNode.getHeight());
			}
		}
		return trHeights;
	}

	public static void getRecipients(Node trTreeSubRoot, List<Node> recipients) {
		if (trTreeSubRoot.isLeaf())
			return;
		Node leftChild = trTreeSubRoot.getChild(0);
		getRecipients(leftChild, recipients);
		if (trTreeSubRoot.getChildCount() > 1 && !trTreeSubRoot.isFake()) {
			Node rightChild = trTreeSubRoot.getChild(1);
			recipients.add(rightChild);
			getRecipients(rightChild, recipients);
		}
		return;
	}

	public static List<Node> getTrueNodeSubTree(Node root) {
		List<Node> subtreeNodes = new ArrayList<>();

		subtreeNodes.add(root);
		for (Node child : Pitchforks.getLogicalChildren(root))
			subtreeNodes.addAll(getTrueNodeSubTree(child));

		return subtreeNodes;
	}

	public static double findMinHeight(Node geneNode, List<Integer> possibleNodeAssignments,
			Integer[] geneNodeAssignment, SpeciesTreeInterface transmissionTree) {
		int parentAssignmentLabel = geneNodeAssignment[geneNode.getNr()];
		if (possibleNodeAssignments.contains(parentAssignmentLabel))
			return geneNode.getHeight();
		Node trNode = transmissionTree.getNode(parentAssignmentLabel);
		while (true) {
			if (possibleNodeAssignments.contains(trNode.getNr()))
				return trNode.getHeight();
			trNode = trNode.getParent();

		}
	}

	public static Boolean fillAndCheck(SpeciesTreeInterface transmissionTree, Node subRoot,
			Integer[] geneTreeNodeAssignment,
			double[] trNodeOccupancy) {

		if (subRoot.isDirty() == Tree.IS_CLEAN
				&& transmissionTree.getNode(geneTreeNodeAssignment[subRoot.getNr()]).isDirty() == Tree.IS_CLEAN) {
			if (!subRoot.isLeaf()) {
				if (!fillAndCheck(transmissionTree, subRoot.getChild(0), geneTreeNodeAssignment,
						trNodeOccupancy))
					return false;
				if (subRoot.getChildCount() > 1
						&& !fillAndCheck(transmissionTree, subRoot.getChild(1), geneTreeNodeAssignment,
								trNodeOccupancy))
					return false;
			}
		} else if (!fillAssignmentAndCheck(transmissionTree, subRoot, geneTreeNodeAssignment,
				trNodeOccupancy))
				return false;


		return true;

	}

	/**
	 * @param transmissionTree
	 * @param subRoot                root at which to start filling the rest of
	 *                               geneTreeNodeAssignment.
	 * @param geneTreeNodeAssignment List of gene tree node assignments in the form
	 *                               geneTreeNodeNr=transmissionTreeNodeNr. Must
	 *                               provide it already filled with leaf node
	 *                               assignments.
	 * @return true if current gene tree is compatible with current transmission
	 *         tree. Fill geneTreeNodeAssignment.
	 */
	public static Boolean fillAssignmentAndCheck(SpeciesTreeInterface transmissionTree, Node subRoot,
			Integer[] geneTreeNodeAssignment,
			double[] trNodeOccupancy) {



		if (!subRoot.isLeaf()) {
			if (!subRoot.isRoot() && greaterWithPrecision(subRoot.getHeight(), subRoot.getParent().getHeight())) {

				throw new RuntimeException("Negative branch length!");
			}

			if (!fillAssignmentAndCheck(transmissionTree, subRoot.getChild(0), geneTreeNodeAssignment,
					trNodeOccupancy))
				return false;
			Node tr1 = transmissionTree
					.getNode(geneTreeNodeAssignment[subRoot.getChild(0).getNr()]);
			while (!tr1.isRoot()) {
				Node tr1ParentNode = tr1.getParent();
				if (greaterOrEqualHeightWithPrecision(tr1ParentNode, subRoot))
					break;
				tr1 = tr1ParentNode;
			}

			if (subRoot.getChildCount() > 1) {
				if (!fillAssignmentAndCheck(transmissionTree, subRoot.getChild(1), geneTreeNodeAssignment,
						trNodeOccupancy))
					return false;
				Node tr2 = transmissionTree
						.getNode(geneTreeNodeAssignment[subRoot.getChild(1).getNr()]);
				while (!tr2.isRoot()) {
					Node tr2ParentNode = tr2.getParent();
					if (greaterOrEqualHeightWithPrecision(tr2ParentNode, subRoot))
						break;
					tr2 = tr2ParentNode;
				}

				if (tr1.getNr() != tr2.getNr())
					return false;
			}
			if (equalHeightWithPrecision(tr1, subRoot) && !subRoot.isLeaf())
				return false;
			boolean recipient = isRecipient(tr1);
			if (!tr1.isRoot() && !recipient && !subRoot.isLeaf()
					&& equalHeightWithPrecision(tr1.getParent(), subRoot))
				return false;
			geneTreeNodeAssignment[subRoot.getNr()] = tr1.getNr();


			// put geneTree lineage lengths per transmission tree
			if (!subRoot.isRoot()) {
				Arrays.fill(trNodeOccupancy, subRoot.getNr() * transmissionTree.getNodeCount(),
						subRoot.getNr() * transmissionTree.getNodeCount() + transmissionTree.getNodeCount(), 0);

				double parentHeight = subRoot.getParent().getHeight();
				
				if (!tr1.isRoot()) {
					if (greaterOrEqualWithPrecision(parentHeight, tr1.getParent().getHeight())) {
						trNodeOccupancy[subRoot.getNr() * transmissionTree.getNodeCount()
								+ tr1.getNr()] += equalWithPrecision(tr1.getParent().getHeight(),
										subRoot.getHeight())
												? 0.0
												: tr1.getParent().getHeight() - subRoot.getHeight();
						Node child = tr1.getParent();
						Node parent = child.getParent();
						while (greaterOrEqualWithPrecision(parentHeight, child.getHeight())) {
							if (!child.isRoot() && greaterOrEqualWithPrecision(parentHeight, parent.getHeight()))
								trNodeOccupancy[subRoot.getNr() * transmissionTree.getNodeCount()
										+ child.getNr()] = equalWithPrecision(parent.getHeight(),
												child.getHeight())
														? 0.0
														: parent.getHeight() - child.getHeight();
							else
								trNodeOccupancy[subRoot.getNr() * transmissionTree.getNodeCount()
										+ child.getNr()] = equalWithPrecision(parentHeight, child.getHeight())
												? 0.0
												: parentHeight - child.getHeight();

							if (child.isRoot())
								break;
							child = child.getParent();
							parent = child.getParent();
						}
					} else {
						trNodeOccupancy[subRoot.getNr() * transmissionTree.getNodeCount()
								+ tr1.getNr()] = equalWithPrecision(parentHeight, subRoot.getHeight())
										? 0.0
										: parentHeight - subRoot.getHeight();
					}
				} else {
					trNodeOccupancy[subRoot.getNr() * transmissionTree.getNodeCount()
							+ tr1.getNr()] = equalWithPrecision(parentHeight, subRoot.getHeight())
									? 0.0
									: parentHeight - subRoot.getHeight();
				}

			}
				

		} else {

			double parentHeight = subRoot.getParent().getHeight();
			Node trLeaf = transmissionTree.getNode(geneTreeNodeAssignment[subRoot.getNr()]);

			Arrays.fill(trNodeOccupancy, subRoot.getNr() * transmissionTree.getNodeCount(),
					subRoot.getNr() * transmissionTree.getNodeCount() + transmissionTree.getNodeCount(), 0);
			if (!trLeaf.isRoot()) {
				if (greaterOrEqualWithPrecision(parentHeight, trLeaf.getParent().getHeight())) {
					trNodeOccupancy[subRoot.getNr() * transmissionTree.getNodeCount()
							+ trLeaf.getNr()] += equalWithPrecision(trLeaf.getParent().getHeight(),
									subRoot.getHeight())
											? 0.0
											: trLeaf.getParent().getHeight() - subRoot.getHeight();
					Node child = trLeaf.getParent();
					Node parent = child.getParent();
					while (greaterOrEqualWithPrecision(parentHeight, child.getHeight())) {
						if (!child.isRoot() && greaterOrEqualWithPrecision(parentHeight, parent.getHeight()))
							trNodeOccupancy[subRoot.getNr() * transmissionTree.getNodeCount()
									+ child.getNr()] = equalWithPrecision(parent.getHeight(), child.getHeight())
											? 0.0
											: parent.getHeight() - child.getHeight();
						else
							trNodeOccupancy[subRoot.getNr() * transmissionTree.getNodeCount()
									+ child.getNr()] = equalWithPrecision(parentHeight, child.getHeight())
											? 0.0
											: parentHeight - child.getHeight();

						if (child.isRoot())
							break;
						child = child.getParent();
						parent = child.getParent();
					}
				} else {
					trNodeOccupancy[subRoot.getNr() * transmissionTree.getNodeCount()
							+ trLeaf.getNr()] = equalWithPrecision(parentHeight, subRoot.getHeight())
									? 0.0
									: parentHeight - subRoot.getHeight();
				}
			} else {
				trNodeOccupancy[subRoot.getNr() * transmissionTree.getNodeCount()
						+ trLeaf.getNr()] = equalWithPrecision(parentHeight, subRoot.getHeight())
								? 0.0
								: parentHeight - subRoot.getHeight();
			}

		}

		try {
			List<Double> trHeights = getTransmissionHeights(transmissionTree);
			if (containsDoubleWithPrecision(trHeights, subRoot.getHeight()) &&
					!equalHeightWithPrecision(subRoot,
							transmissionTree.getNode(geneTreeNodeAssignment[subRoot.getNr()]).getParent()))
			return false;
		} catch (Exception e) {
			throw new RuntimeException("Exception while checking the gene tree");
		}
		return true;
	}

	public static boolean isRecipient(Node trTreeNode) {
		return (!trTreeNode.isRoot()
				&& trTreeNode.getParent().getChild(0) != trTreeNode && !trTreeNode.getParent().isFake());
	}

	public static boolean equalHeightWithPrecision(Node n1, Node n2) {
		if (Math.abs(n1.getHeight() - n2.getHeight()) <= globalPrecisionThreshold)
			return true;
		return false;
	}

	public static boolean equalWithPrecision(Double d1, Double d2) {
		if (Math.abs(d1 - d2) <= globalPrecisionThreshold)
			return true;
		return false;
	}

	public static boolean greaterOrEqualHeightWithPrecision(Node n1, Node n2) {
		if (n1.getHeight() > n2.getHeight() + globalPrecisionThreshold || equalHeightWithPrecision(n1, n2))
			return true;
		return false;
	}

	public static boolean greaterHeightNode(Node n1, Node n2) {
		if (n1.getHeight() > n2.getHeight() + globalPrecisionThreshold)
			return true;
		return false;
	}

	public static boolean greaterOrEqualWithPrecision(Double d1, Double d2) {
		if (greaterWithPrecision(d1, d2) || equalWithPrecision(d1, d2))
			return true;
		return false;
	}

	public static boolean containsDoubleWithPrecision(Set<Double> set, Double d) {
		for (Double l : set) {
			if (equalWithPrecision(d, l))
				return true;
		}

		return false;
	}

	public static boolean greaterWithPrecision(Double d1, Double d2) {
		if (d1 > d2 + globalPrecisionThreshold)
			return true;
		return false;
	}

	public static boolean containsDoubleWithPrecision(List<Double> list, Double d) {
		for (int i=0; i<list.size(); i++) {
			if (equalWithPrecision(list.get(i), d))
				return true;
		}
		return false;
	}

	public static double round(double value, int places) {
		if (places < 0)
			throw new IllegalArgumentException();

		BigDecimal bd = BigDecimal.valueOf(value);
		bd = bd.setScale(places, RoundingMode.HALF_UP);
		return bd.doubleValue();
	}

	public static void exchangeNodesRandomDirection(Node i, Node j,
			Node p, Node jP) {
		// precondition p -> i & jP -> j
		replaceRandomDirection(p, i, j);
		replaceRandomDirection(jP, j, i);
		// postcondition p -> j & p -> i
	}

	public static void replaceRandomDirection(final Node node, final Node child, final Node replacement) {
		Node otherChild = getOtherChild(node, child);
		node.removeChild(otherChild);
		node.removeChild(child);
		if (Randomizer.nextBoolean()) {
			node.addChild(replacement);
			node.addChild(otherChild);
		} else {
			node.addChild(otherChild);
			node.addChild(replacement);
		}
		node.makeDirty(Tree.IS_FILTHY);
		replacement.makeDirty(Tree.IS_FILTHY);
	}

	public static void exchangeNodesKeepDirection(Node i, Node j,
			Node p, Node jP) {
		// precondition p -> i & jP -> j
		replaceNodeKeepDirection(p, i, j);
		replaceNodeKeepDirection(jP, j, i);
		// postcondition p -> j & p -> i
	}

	public static void replaceNodeKeepDirection(final Node node, final Node child, final Node replacement) {
			boolean left = node.getLeft().getNr() == child.getNr();
			Node otherChild = getOtherChild(node, child);
			boolean fake = otherChild.isDirectAncestor();
			node.removeChild(otherChild);
			node.removeChild(child);
			if (left || fake) {
				node.addChild(replacement);
				node.addChild(otherChild);
			} else {
				node.addChild(otherChild);
				node.addChild(replacement);
			}
			node.makeDirty(Tree.IS_FILTHY);
			replacement.makeDirty(Tree.IS_FILTHY);
		}

		/**
		 * Divide the branch which is NOT a root branch with a single child node.
		 * 
		 * @param node
		 * @param child
		 * @param dividerNode
		 * @return
		 */
		public static void divideBranch(final Node node, final Node child, final Node dividerNode) {
			if (node.getChildCount() > 1) {
				replaceNodeKeepDirection(node, child, dividerNode);
			} else {
				node.removeChild(child);
				node.addChild(dividerNode);
				node.makeDirty(Tree.IS_FILTHY);
				dividerNode.makeDirty(Tree.IS_FILTHY);
			}
			dividerNode.addChild(child);
			child.makeDirty(Tree.IS_FILTHY);
		}

	/**
	 * @param parent the parent
	 * @param child  the child that you want the sister of
	 * @return the other child of the given parent.
	 */
	protected static Node getOtherChild(final Node parent, final Node child) {
		if (parent.getLeft().getNr() == child.getNr()) {
			return parent.getRight();
		} else {
			return parent.getLeft();
		}
	}

	/**
	 * COPIED from BDMM prime package by T. G. Vaughan on 2022-03-12
	 * 
	 * Apply node numbers to internal nodes below and including subtreeRoot. Numbers
	 * are applied postorder, so parents always have larger numbers than their
	 * children and the root has the hightest number.
	 *
	 * @param subtreeRoot root of subtree
	 * @param nextNumber  next number to be used
	 * @return next number to be used on another part of the tree.
	 */
	public static int numberInternalNodesOnSubtree(Node subtreeRoot, int nextNumber) {

		if (subtreeRoot.isLeaf())
			return nextNumber;

		for (Node child : subtreeRoot.getChildren())
			nextNumber = numberInternalNodesOnSubtree(child, nextNumber);

		subtreeRoot.setNr(nextNumber);

		return nextNumber + 1;
	}

	public static void getGroupAndLogicalNotHiddenChildren(Node node, List<Node> group, List<Node> logicalChildren) {
		for (Node child : node.getChildren()) {
			if (! (child.getChildCount()==2 && child.isFake())){
				if (child.getHeight() == node.getHeight()) {
					if (group != null)
						group.add(child);
					getGroupAndLogicalNotHiddenChildren(child, group, logicalChildren);
				} else {
					if (logicalChildren != null)
						logicalChildren.add(child);
				}
			}
		}
	}

	/**
	 * Get logical children descending from logical node groupRoot.
	 *
	 * @param groupRoot logical node in Tree
	 * @return newly created list of child nodes.
	 */
	public static List<Node> getLogicalNotHiddenChildren(Node groupRoot) {
		List<Node> logicalChildren = new ArrayList<>();

		getGroupAndLogicalNotHiddenChildren(groupRoot, null, logicalChildren);

		return logicalChildren;
	}

}
