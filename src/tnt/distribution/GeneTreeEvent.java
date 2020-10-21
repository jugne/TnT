package tnt.distribution;

import java.util.ArrayList;
import java.util.List;

import beast.evolution.tree.Node;

public class GeneTreeEvent {
	public enum GeneTreeEventType {
		SAMPLE, BIFURCATION, MULTIFURCATION, MOCK
	}

	public GeneTreeEventType type;
	public double time;

	// Only used for multifurcation events,
	// since they are encoded as sequence of bifurcations
	public int fakeBifCount = 0;

	// Only used for multifurcation events,
	// how many separate bif/multi-furcating events happen during bottleneck
	public int multiCoalCount = 0;

	public List<Integer> multiCoalSize = new ArrayList<>();

	public int lineages;

	public Node node;

	// numbers of tree nodes involved in the event
	public List<Integer> nodesInEventNr;
}
