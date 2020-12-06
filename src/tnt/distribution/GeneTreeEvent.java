package tnt.distribution;

import java.util.ArrayList;
import java.util.List;

import beast.evolution.tree.Node;


public class GeneTreeEvent {
	// bifurcation events are only those, that have SINGLE bifurcation happening on
	// a tree at this event time point
	// multifurcation event includes also situations with MORE THAN ONE bifurcation
	// happening at the exactly same time (multiple mergers)
	public enum GeneTreeEventType {
		SAMPLE, BIFURCATION, MULTIFURCATION, MOCK
	}

	public GeneTreeEventType type;
	public double time;

	// Only used for MULTIFURCATION events,
	// since they are encoded as sequence of bifurcations.
	// Fake bifurcation is a bifurcating node that has a parent edge with length 0.
	public int fakeBifCount = 0;

	// Only used for MULTIFURCATION events,
	// how many separate bif/multi-furcating events happen during bottleneck
	public int multiCoalCount = 0;

	public List<Integer> multiCoalSize = new ArrayList<>();

	public int lineages;

	public Node node;

	// numbers of tree nodes involved in the event
	public List<Integer> nodesInEventNr;
}
