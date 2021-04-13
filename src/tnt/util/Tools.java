package tnt.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import pitchfork.Pitchforks;
import starbeast2.SpeciesTreeInterface;

public class Tools {

	public static boolean isMultiMerger(List<Node> allLogicalNodes, Node logicalNode) {
		for (Node n : allLogicalNodes) {
			if (n.getHeight() == logicalNode.getHeight() && n.getNr() != logicalNode.getNr())
				return true;
		}
		return false;
	}

	public static List<Node> getGeneNodesAtTransmission(List<Node> nodeList, List<Double> transmissionHeights) {
		List<Node> tmp = new ArrayList<>();
		for (Node n : nodeList) {
			if (transmissionHeights.contains(n.getHeight()))
				tmp.add(n);
		}

		return tmp;
	}


	public static List<Node> getGeneNodesNotAtTransmission(List<Node> nodeList, List<Double> transmissionHeights) {
		List<Node> tmp = new ArrayList<>();
		for (Node n : nodeList) {
			if (!transmissionHeights.contains(n.getHeight()))
				tmp.add(n);
		}

		return tmp;
	}

	public static List<Node> getGeneNodesWithParentNotAtTransmission(List<Node> nodeList,
			List<Double> transmissionHeights) {
		List<Node> tmp = new ArrayList<>();
		for (Node n : nodeList) {
			if (n.isRoot() || !transmissionHeights.contains(n.getParent().getHeight()))
				tmp.add(n);
		}

		return tmp;
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
			else if (transmissionHeights.contains(n.getParent().getHeight()))
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
			if (!subRoot.isRoot() && subRoot.getParent().getHeight() < subRoot.getHeight()) {

				throw new RuntimeException("Negative branch length!");
			}

			if (!fillAssignmentAndCheck(transmissionTree, subRoot.getChild(0), geneTreeNodeAssignment,
					trNodeOccupancy))
				return false;
			Node tr1 = transmissionTree
					.getNode(geneTreeNodeAssignment[subRoot.getChild(0).getNr()]);
			while (!tr1.isRoot()) {
				Node tr1ParentNode = tr1.getParent();
				if (tr1ParentNode.getHeight() >= subRoot.getHeight())
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
					if (tr2ParentNode.getHeight() >= subRoot.getHeight())
						break;
					tr2 = tr2ParentNode;
				}

				if (tr1.getNr() != tr2.getNr())
					return false;
			}
			if (tr1.getHeight() == subRoot.getHeight() && !subRoot.isLeaf())
				return false;
			boolean recipient = isRecipient(tr1);
			if (!tr1.isRoot() && !recipient && !subRoot.isLeaf() && tr1.getParent().getHeight() == subRoot.getHeight())
				return false;
			geneTreeNodeAssignment[subRoot.getNr()] = tr1.getNr();


			// put geneTree lineage lengths per transmission tree
			if (!subRoot.isRoot()) {
				Arrays.fill(trNodeOccupancy, subRoot.getNr() * transmissionTree.getNodeCount(),
						subRoot.getNr() * transmissionTree.getNodeCount() + transmissionTree.getNodeCount(), 0);

				double parentHeight = subRoot.getParent().getHeight();
				
				if (!tr1.isRoot()) {
					if (parentHeight >= tr1.getParent().getHeight()) {
						trNodeOccupancy[subRoot.getNr() * transmissionTree.getNodeCount() + tr1.getNr()] += tr1
								.getParent().getHeight() - subRoot.getHeight();
						Node child = tr1.getParent();
						Node parent = child.getParent();
						while (parentHeight >= child.getHeight()) {
							if (!child.isRoot() && parentHeight >= parent.getHeight())
								trNodeOccupancy[subRoot.getNr() * transmissionTree.getNodeCount()
										+ child.getNr()] = parent.getHeight() - child.getHeight();
							else
								trNodeOccupancy[subRoot.getNr() * transmissionTree.getNodeCount()
										+ child.getNr()] = parentHeight - child.getHeight();

							if (child.isRoot())
								break;
							child = child.getParent();
							parent = child.getParent();
						}
					} else {
						trNodeOccupancy[subRoot.getNr() * transmissionTree.getNodeCount()
								+ tr1.getNr()] = parentHeight - subRoot.getHeight();
					}
				} else {
					trNodeOccupancy[subRoot.getNr() * transmissionTree.getNodeCount()
							+ tr1.getNr()] = parentHeight - subRoot.getHeight();
				}

			}
				

		} else {

			double parentHeight = subRoot.getParent().getHeight();
			Node trLeaf = transmissionTree.getNode(geneTreeNodeAssignment[subRoot.getNr()]);

			Arrays.fill(trNodeOccupancy, subRoot.getNr() * transmissionTree.getNodeCount(),
					subRoot.getNr() * transmissionTree.getNodeCount() + transmissionTree.getNodeCount(), 0);
			if (!trLeaf.isRoot()) {
				if (parentHeight >= trLeaf.getParent().getHeight()) {
					trNodeOccupancy[subRoot.getNr() * transmissionTree.getNodeCount() + trLeaf.getNr()] += trLeaf
							.getParent().getHeight() - subRoot.getHeight();
					Node child = trLeaf.getParent();
					Node parent = child.getParent();
					while (parentHeight >= child.getHeight()) {
						if (!child.isRoot() && parentHeight >= parent.getHeight())
							trNodeOccupancy[subRoot.getNr() * transmissionTree.getNodeCount()
									+ child.getNr()] = parent.getHeight() - child.getHeight();
						else
							trNodeOccupancy[subRoot.getNr() * transmissionTree.getNodeCount()
									+ child.getNr()] = parentHeight - child.getHeight();

						if (child.isRoot())
							break;
						child = child.getParent();
						parent = child.getParent();
					}
				} else {
					trNodeOccupancy[subRoot.getNr() * transmissionTree.getNodeCount()
							+ trLeaf.getNr()] = parentHeight - subRoot.getHeight();
				}
			} else {
				trNodeOccupancy[subRoot.getNr() * transmissionTree.getNodeCount()
						+ trLeaf.getNr()] = parentHeight - subRoot.getHeight();
			}

		}

		try {
			List<Double> trHeights = getTransmissionHeights(transmissionTree);
			if (trHeights.contains(subRoot.getHeight()) &&
					subRoot.getHeight() != transmissionTree.getNode(geneTreeNodeAssignment[subRoot.getNr()]).getParent()
						.getHeight())
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
}
