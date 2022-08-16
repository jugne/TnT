package tnt.operators;

import java.util.ArrayList;
import java.util.List;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import tnt.distribution.GeneTreeIntervals;
import tnt.util.Tools;

/**
 *@author Alexandra Gavryushkina
 */
public class SAWilsonBalding extends TreeOperator {

    public Input<RealParameter> rInput =
            new Input<RealParameter>("removalProbability", "The probability of an individual to become noninfectious immediately after the sampling");
	public Input<List<GeneTreeIntervals>> geneTreeIntervalsInput = new Input<>("geneTreeIntervals",
			"intervals for a gene tree", new ArrayList<>(), Validate.REQUIRED);

	public Input<Boolean> allowSACreationInput = new Input<Boolean>("allowSaCreation",
			"If the operator should be able to create SA nodes. Default 'true'. "
					+ "Only applicable if removalProbability >0.",
			true);

    @Override
    public void initAndValidate() {
    }

    /**
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {

        Tree tree = treeInput.get(this);

		int oldSACount = tree.getDirectAncestorNodeCount();

        //double x0 = 10;

        double oldMinAge, newMinAge, newRange, oldRange, newAge, fHastingsRatio, DimensionCoefficient;
        int newDimension, oldDimension;

		// get fit nodes as list
		Integer[] fitNodesNrs = Tools.getTrNodeNrsNotTransmissionOnGenes(tree, geneTreeIntervalsInput.get(), true);
        
//		Integer[] fitNodesNrs2 = getTrNodeNrsNotTransmissionOnGenesAfter(tree);

        // choose a random node avoiding root and leaves that are direct ancestors
		int nodeCount = tree.getNodeCount();
		int nodeCountFrom = fitNodesNrs.length;
//		int nodeCountFrom2 = fitNodesNrs2.length;
		if (nodeCountFrom == 1)
			return Double.NEGATIVE_INFINITY;
        Node i;

        do {
//            i = tree.getNode(Randomizer.nextInt(nodeCount));
			int nr = fitNodesNrs[Randomizer.nextInt(nodeCountFrom)];
			i = tree.getNode(nr);
		} while (i.isRoot() || i.isDirectAncestor());

//		if (i.isDirectAncestor())
//			throw new RuntimeException("Direct ancestor node selected in wilson balding");

        Node iP = i.getParent();
        Node CiP;
        if (iP.getLeft().getNr() == i.getNr()) {
            CiP = iP.getRight();
        } else {
            CiP = iP.getLeft();
        }

        // make sure that there is at least one candidate edge to attach node iP to
		if (iP.getParent() == null && Tools.greaterOrEqualHeightWithPrecision(i, CiP)) {
            return Double.NEGATIVE_INFINITY;
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
//				j = tree.getNode(nodeNumber);
				j = tree.getExternalNodes().get(nodeNumber - nodeCount);
                jP = j.getParent();
                newParentHeight = j.getHeight();
                attachingToLeaf = true;
                //adjacentLeaf = (iP.getNr() == j.getNr());
            }
		} while (j.isDirectAncestor() || (Tools.greaterOrEqualWithPrecision(i.getHeight(), newParentHeight))
				|| (i.getNr() == j.getNr()) || adjacentEdge /* || adjacentLeaf */);


        if (attachingToLeaf && iP.getNr() == j.getNr()) {
            System.out.println("Proposal failed because j = iP");
            return Double.NEGATIVE_INFINITY;
        }

        if (jP != null && jP.getNr() == i.getNr()) {
            System.out.println("Proposal failed because jP = i. Heights of i = " + i.getHeight() + " Height of jP = " + jP.getHeight());
            return Double.NEGATIVE_INFINITY;
        }

//		if (nodeCountFrom2 != nodeCountFrom)
//			System.out.println();
		oldDimension = nodeCountFrom - tree.getDirectAncestorNodeCount() - 1;

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
			if (right && !j.isLeaf()) {
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

		int newSACount = tree.getDirectAncestorNodeCount();
		if (!allowSACreationInput.get() && oldSACount != newSACount) {
			return Double.NEGATIVE_INFINITY;
		}

//		newDimension = nodeCountFrom - tree.getDirectAncestorNodeCount() - 1;
//		if (nodeCountFrom != getTrNodeNrsNotTransmissionOnGenesAfter(tree).length)
//			System.out.println();
		newDimension = nodeCountFrom - tree.getDirectAncestorNodeCount() - 1;
        DimensionCoefficient = (double) oldDimension / newDimension;

        fHastingsRatio = Math.abs(DimensionCoefficient * newRange / oldRange);

//		System.out.println(PiP.getNr());
//		System.out.println(i.getNr());
        return Math.log(fHastingsRatio);

    }

//	private Integer[] getTrNodeNrsNotTransmissionOnGenes(Tree transmissionTree) {
//		Set<Integer> fitNodeNrs = new HashSet<Integer>();
//		final List<GeneTreeIntervals> intervalsLis = geneTreeIntervalsInput.get();
//
//		for (Node n : transmissionTree.getNodesAsArray()) {
//			if (n.isRoot()) {
//				fitNodeNrs.add(n.getNr());
//				continue;
//			}
//
//			int recipientNr = n.getParent().getChild(1).getNr();
//			boolean addToFit = true;
//
//			for (GeneTreeIntervals intervals : intervalsLis) {
//				HashMap<Integer, List<GeneTreeEvent>> eventList = intervals.getGeneTreeEventList();
//				List<GeneTreeEvent> eventsPerTrNode = eventList.get(recipientNr);
//				GeneTreeEvent lastEvent = eventsPerTrNode.get(eventsPerTrNode.size() - 1);
//
//				if (!n.getParent().isFake()
//						&& Tools.equalWithPrecisionDouble(n.getParent().getHeight(), lastEvent.time)) {
//					addToFit = false;
//					break;
//				}
//			}
//
//			if (addToFit) {
//				fitNodeNrs.add(recipientNr);
//				fitNodeNrs.add(n.getParent().getChild(0).getNr());
//			}
//		}
//
//		return fitNodeNrs.toArray(new Integer[0]);
//	}

//	private Integer[] getTrNodeNrsNotTransmissionOnGenesAfter(Tree transmissionTree) {
//		Set<Double> unfitHeights = new HashSet<Double>();
//		final List<Double> trHeights = Tools.getTransmissionHeights((TransmissionTree) transmissionTree);
//		Double trRootHeight = transmissionTree.getRoot().getHeight();
//		final List<GeneTreeIntervals> intervalsLis = geneTreeIntervalsInput.get();
//
//		for (GeneTreeIntervals intervals : intervalsLis) {
//			List<Node> nodes = Pitchforks.getTrueInternalNodes(intervals.geneTreeInput.get());
//			for (Node n : nodes) {
//				if (!n.isLeaf() && trHeights.contains(n.getHeight())) {
//					unfitHeights.add(n.getHeight());
//				}
//			}
//		}
//
////		unfitHeights.remove(trRootHeight);
//		Set<Integer> fitNodeNrs = new HashSet<Integer>();
//		for (Node trNode : transmissionTree.getNodesAsArray()) {
//			if (trNode.isRoot())
//				fitNodeNrs.add(trNode.getNr());
//			else if (trNode.getParent().isFake() || !unfitHeights.contains(trNode.getParent().getHeight()))
//				fitNodeNrs.add(trNode.getNr());
//		}
//
//		return fitNodeNrs.toArray(new Integer[0]);
//	}


}
