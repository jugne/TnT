package tnt.util;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import beast.evolution.tree.Node;
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

	public static List<Node> getGeneNodesWithParentsAtTransmission(List<Node> nodeList,
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
			HashMap<Integer, Integer> geneTreeNodeAssignment) {
		if (!subRoot.isLeaf()) {
			if (!fillAssignmentAndCheck(transmissionTree, subRoot.getChild(0), geneTreeNodeAssignment))
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
				if (!fillAssignmentAndCheck(transmissionTree, subRoot.getChild(1), geneTreeNodeAssignment))
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
			boolean recipient = isRecipient(tr1);
			if (!tr1.isRoot() && !recipient && !subRoot.isLeaf() && tr1.getParent().getHeight() == subRoot.getHeight())
				return false;
			geneTreeNodeAssignment.put(subRoot.getNr(), tr1.getNr());
		}
		return true;
	}

	public static boolean isRecipient(Node trTreeNode) {
		return (!trTreeNode.isRoot()
				&& trTreeNode.getParent().getChild(0) != trTreeNode && !trTreeNode.getParent().isFake());
	}
}
