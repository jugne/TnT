package tnt.likelihood;

import beast.evolution.tree.Node;

public class GeneTreeEvent {
	public enum GeneTreeEventType {
		SAMPLE, BIFURCATION, MULTIFURCATION
	}

	public GeneTreeEventType type;
	public double time;

	// Only used for multifurcation events,
	// since they are encoded as sequence of bifurcations
	public int fakeBifCount = 0;

	// Only used for multifurcation events,
	// how many separate bif/multi-furcating events happen during bottleneck
	public int multiCoalCount = 0;

	public int lineages;

	public Node node;
}
