package tnt.util;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import beast.evolution.tree.Node;

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

	public static <T> boolean listEqualsIgnoreOrder(List<T> list1, List<T> list2) {
		return new HashSet<>(list1).equals(new HashSet<>(list2));
	}
}
