package tnt.operators;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import tnt.distribution.GeneTreeEvent;
import tnt.distribution.GeneTreeIntervals;

/**
 *@author Alexandra Gavryushkina
 */
public class SAWilsonBalding extends TreeOperator {

    public Input<RealParameter> rInput =
            new Input<RealParameter>("removalProbability", "The probability of an individual to become noninfectious immediately after the sampling");
	public Input<List<GeneTreeIntervals>> geneTreeIntervalsInput = new Input<>("geneTreeIntervals",
			"intervals for a gene tree", new ArrayList<>(), Validate.REQUIRED);

    @Override
    public void initAndValidate() {
    }

    /**
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {

        Tree tree = treeInput.get(this);

        //double x0 = 10;

		double oldMinAge, newMinAge, newRange, oldRange, newAge, fHastingsRatio, DimensionCoefficient,
				eventTimeCoefficient;
        int newDimension, oldDimension;
        
        // weather to choose node to detach (i) which parent (iP) coinsides with gene tree event
        boolean detachFromGeneTreeEventTime = Randomizer.nextBoolean();
        // weather to choose node to attach (j) which parent (jP) coinsides with gene tree event
        boolean attachAtGeneTreeEventTime = Randomizer.nextBoolean();

		// get nodes that do not have parents that coincide with transmission induced
		// events on gene trees
		Set<Integer> nodeNrsWithParentNOTAtGeneTreeTransmission = getTrNodesWithParentNotAtGeneNodeTimes(tree);
        
        // choose a random node avoiding root and leaves that are direct ancestors
		int nodeCount = tree.getNodeCount();
		int nodeCountFrom;
		if (!detachFromGeneTreeEventTime)
			nodeCountFrom = nodeNrsWithParentNOTAtGeneTreeTransmission.size();
		else
			nodeCountFrom = nodeCount - nodeNrsWithParentNOTAtGeneTreeTransmission.size();
		
		if (nodeCountFrom == 1)
			return Double.NEGATIVE_INFINITY;
        
		
		Node i;
        do {
            i = tree.getNode(Randomizer.nextInt(nodeCount));
//			int nr = nodeNrsWithParentAtGeneTreeTransmission[Randomizer.nextInt(nodeCountFrom)];
//			i = tree.getNode(nr);
		} while (((detachFromGeneTreeEventTime && nodeNrsWithParentNOTAtGeneTreeTransmission.contains(i.getNr()))
				|| (!detachFromGeneTreeEventTime && !nodeNrsWithParentNOTAtGeneTreeTransmission.contains(i.getNr())))
				&& (i.isRoot() || i.isDirectAncestor()));


        Node iP = i.getParent();
        Node CiP;
        if (iP.getLeft().getNr() == i.getNr()) {
            CiP = iP.getRight();
        } else {
            CiP = iP.getLeft();
        }

        // make sure that there is at least one candidate edge to attach node iP to
        if (iP.getParent() == null && CiP.getHeight() <= i.getHeight()) {
            return Double.NEGATIVE_INFINITY;
        }

		// check for all gene trees
		final List<GeneTreeIntervals> intervalsList = geneTreeIntervalsInput.get();
		// ways to choose attachment time if detachFromGeneTreeEventTime
		int oldPossibleAttachEventTimesCount = 0;
		if (detachFromGeneTreeEventTime) {

			// events, excluding sample and mock, corresponding to i, iP and CiP
			// transmission tree nodes in all gene trees
			List<GeneTreeEvent> oldEventsToAttach = null;
			Set<Double> oldPossibleAttachEventTimes = new HashSet<Double>();
			for (GeneTreeIntervals intervals : intervalsList) {
				HashMap<Integer, List<GeneTreeEvent>> eventList = intervals.getGeneTreeEventList();
				oldEventsToAttach = eventList.get(i.getNr());
				oldEventsToAttach.addAll(eventList.get(iP.getNr()));
				oldEventsToAttach.addAll(eventList.get(CiP.getNr()));

//				oldEventsToAttach.removeIf(e -> e.type == GeneTreeEventType.SAMPLE || e.type == GeneTreeEventType.MOCK);
			}
			for (GeneTreeEvent e : oldEventsToAttach) {
				if (oldPossibleAttachEventTimes.contains(e.time))
					continue;
				oldPossibleAttachEventTimes.add(e.time);
			}
			oldPossibleAttachEventTimesCount = oldPossibleAttachEventTimes.size();
		}

        // choose another random node to insert i above or to attach i to this node if it is a leaf
        Node j;
        Node jP;

        final int leafNodeCount = tree.getLeafNodeCount();

        if (leafNodeCount != tree.getExternalNodes().size()) {
            System.out.println("node counts are incorrect. NodeCount = " + nodeCount + " leafNodeCount = " + leafNodeCount + " exteranl node count = " + tree.getExternalNodes().size());
        }

        // make sure that the target branch <jP, j> or target leaf j is above the subtree being moved

        int nodeNumber;
        double newParentHeight;
        boolean attachingToLeaf;
        boolean adjacentEdge;
        //boolean adjacentLeaf;
        do {
            adjacentEdge = false;
            //adjacentLeaf = false;
            nodeNumber = Randomizer.nextInt(nodeCount + leafNodeCount);
            if (nodeNumber < nodeCount) {
                j = tree.getNode(nodeNumber);
                jP = j.getParent();
                if (jP != null)
                    newParentHeight = jP.getHeight();
                else newParentHeight = Double.POSITIVE_INFINITY;
                if (!CiP.isDirectAncestor())
                    adjacentEdge = (CiP.getNr() == j.getNr() || iP.getNr() == j.getNr());
                attachingToLeaf = false;
            } else {
                j = tree.getExternalNodes().get(nodeNumber - nodeCount);
                jP = j.getParent();
                newParentHeight = j.getHeight();
                attachingToLeaf = true;
                //adjacentLeaf = (iP.getNr() == j.getNr());
            }
        } while (j.isDirectAncestor() || (newParentHeight <= i.getHeight()) || (i.getNr() == j.getNr()) || adjacentEdge /*|| adjacentLeaf */);


        if (attachingToLeaf && iP.getNr() == j.getNr()) {
            System.out.println("Proposal failed because j = iP");
            return Double.NEGATIVE_INFINITY;
        }

        if (jP != null && jP.getNr() == i.getNr()) {
            System.out.println("Proposal failed because jP = i. Heights of i = " + i.getHeight() + " Height of jP = " + jP.getHeight());
            return Double.NEGATIVE_INFINITY;
        }

		// total transmission tree node count - number of direct ancestors - root
		oldDimension = nodeCountFrom - tree.getDirectAncestorNodeCount() - 1;

		if (attachAtGeneTreeEventTime) {
			// events, excluding sample and mock, corresponding to i, iP and CiP
			// transmission tree nodes in all gene trees
			List<GeneTreeEvent> newEventsToAttach = null;
			Set<Double> oldPossibleAttachEventTimes = new HashSet<Double>();
			for (GeneTreeIntervals intervals : intervalsList) {
				HashMap<Integer, List<GeneTreeEvent>> eventList = intervals.getGeneTreeEventList();
				newEventsToAttach = eventList.get(i.getNr());
				oldEventsToAttach.addAll(eventList.get(iP.getNr()));
				oldEventsToAttach.addAll(eventList.get(CiP.getNr()));

//				oldEventsToAttach.removeIf(e -> e.type == GeneTreeEventType.SAMPLE || e.type == GeneTreeEventType.MOCK);
			}
			for (GeneTreeEvent e : oldEventsToAttach) {
				if (oldPossibleAttachEventTimes.contains(e.time))
					continue;
				oldPossibleAttachEventTimes.add(e.time);
			}
			oldPossibleAttachEventTimesCount = oldPossibleAttachEventTimes.size();
		}

        //Hastings numerator calculation + newAge of iP
        if (attachingToLeaf) {
            newRange = 1;
            newAge = j.getHeight();
        } else {
            if (jP != null) {
                newMinAge = Math.max(i.getHeight(), j.getHeight());
                newRange = jP.getHeight() - newMinAge;
                newAge = newMinAge + (Randomizer.nextDouble() * newRange);
            } else {
                double randomNumberFromExponential;
                randomNumberFromExponential = Randomizer.nextExponential(1);
                //newRange = x0 - j.getHeight();
                //randomNumberFromExponential = Randomizer.nextDouble() * newRange;
                newRange = Math.exp(randomNumberFromExponential);
                newAge = j.getHeight() + randomNumberFromExponential;
            }
        }

        Node PiP = iP.getParent();

        //Hastings denominator calculation
        if (CiP.isDirectAncestor()) {
            oldRange = 1;
        }
        else {
            oldMinAge = Math.max(i.getHeight(), CiP.getHeight());
            if (PiP != null) {
                oldRange = PiP.getHeight() - oldMinAge;
            } else {
                oldRange = Math.exp(iP.getHeight() - oldMinAge);
                //oldRange = x0 - oldMinAge;
            }
        }

        //update
        if (iP.getNr() != j.getNr() && CiP.getNr() != j.getNr()) {
			Node otherIPChild = getOtherChild(iP, CiP);
			boolean right = iP.getRight().getNr() == otherIPChild.getNr();
            iP.removeChild(CiP); //remove <iP, CiP>


            if (PiP != null) {
				////
				boolean left = PiP.getLeft().getNr() == iP.getNr();
				Node otherChild = getOtherChild(PiP, iP);
				PiP.removeChild(otherChild);
				PiP.removeChild(iP);
				if (left) {
					PiP.addChild(CiP);
					PiP.addChild(otherChild);
				} else {
					PiP.addChild(otherChild);
					PiP.addChild(CiP);
				}

//                PiP.removeChild(iP);   // remove <PiP,iP>
//                PiP.addChild(CiP);   // add <PiP, CiP>
                PiP.makeDirty(Tree.IS_FILTHY);
                CiP.makeDirty(Tree.IS_FILTHY);
            } else {
                CiP.setParent(null); // completely remove <iP, CiP>
                tree.setRootOnly(CiP);
            }

            if (jP != null) {
				replace(jP, j, iP);

//                jP.removeChild(j);  // remove <jP, j>
//                jP.addChild(iP);   // add <jP, iP>
                jP.makeDirty(Tree.IS_FILTHY);
            } else {
                iP.setParent(null); // completely remove <PiP, iP>
                tree.setRootOnly(iP);
            }
			iP.removeChild(otherIPChild);
			if (right) {
				iP.addChild(j);
				iP.addChild(otherIPChild);
			} else {
				iP.addChild(otherIPChild);
				iP.addChild(j);
			}

            iP.makeDirty(Tree.IS_FILTHY);
            j.makeDirty(Tree.IS_FILTHY);
        }
        iP.setHeight(newAge);

        //make sure that either there are no direct ancestors or r<1
        if ((rInput.get() != null) && (tree.getDirectAncestorNodeCount() > 0 && rInput.get().getValue() == 1))  {
            return Double.NEGATIVE_INFINITY;
        }

		newDimension = nodeCountFrom - tree.getDirectAncestorNodeCount() - 1;
        DimensionCoefficient = (double) oldDimension / newDimension;

        fHastingsRatio = Math.abs(DimensionCoefficient * newRange / oldRange);

        return Math.log(fHastingsRatio);

    }

	/**
	 * Gets the Transmission Tree node numbers for nodes that have parent nodes NOT
	 * coinciding with transmission induced events on Gene Trees.
	 *
	 * @param transmissionTree the transmission tree
	 * @return the transmission node numbers with parent nodes not at transmission
	 *         induced gene node times
	 */
	private Set<Integer> getTrNodesWithParentNotAtGeneNodeTimes(Tree transmissionTree) {
		Set<Integer> fitNodeNrs = new HashSet<Integer>();

		// check for all gene trees
		final List<GeneTreeIntervals> intervalsLis = geneTreeIntervalsInput.get();

		// root does not have a parent node
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
				GeneTreeEvent lastEvent = eventsPerTrNode.get(eventsPerTrNode.size() - 1);

				// fake event coinside with gene tree events, but because if sampling, not
				// transmission
				if (!n.getParent().isFake() && n.getParent().getHeight() == lastEvent.time) {
					addToFit = false;
					break;
				}
			}

			if (addToFit) {
				fitNodeNrs.add(recipientNr);
				fitNodeNrs.add(n.getParent().getChild(0).getNr());
			}
		}
		return fitNodeNrs;
	}


}
