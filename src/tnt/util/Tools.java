package tnt.util;

import java.util.ArrayList;
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
}
