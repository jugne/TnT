package tnt.operators;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Original Exchange operator extended:
 * For trees with sampled ancestors by Alexandra Gavryushkina.
 * For transmission trees by Ugne Stolz.
 */


import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import tnt.distribution.GeneTreeIntervals;
import tnt.util.Tools;

@Description("Implement Narrow and Wide Exchange for sampled ancestor trees." +
        "Narrow move chooses a random internal node (not a fake node) with two non-leaf children." +
        "Then it takes the older child of this node and exchange one of its children (or just a child" +
        "if there is only one) with the younger child. Wide remains the same as for regular trees.")
public class TnTExchange extends TreeOperator {

//	final public Input<Tree> treeInput = new Input<>("tree",
//			"beast.tree on which this operation is performed",
//			Validate.REQUIRED);

	final public Input<Boolean> isNarrowInput = new Input<>("isNarrow",
			"if true (default) a narrow exchange is performed, otherwise a wide exchange", true);

	public Input<List<GeneTreeIntervals>> geneTreeIntervalsInput = new Input<>("geneTreeIntervals",
			"intervals for a gene tree", new ArrayList<>(), Validate.REQUIRED);

	@Override
	public void initAndValidate() {
	}

	/**
	 * override this for proposals,
	 *
	 * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should
	 *         not be accepted *
	 */
	@Override
	public double proposal() {
		final Tree tree = treeInput.get(this);

		double logHastingsRatio = 0;

		if (isNarrowInput.get()) {
			logHastingsRatio = narrow(tree);
		} else {
			logHastingsRatio = wide(tree);
		}

		return logHastingsRatio;
	}

	/**
	 * @param tree
	 * @return log of hastings ratio for the narrow move
	 */
	public double narrow(final Tree tree) {

        final int nodeCount = tree.getNodeCount();

        //make sure that there are at least two distinct non-root nodes which are not direct ancestors.
        if (nodeCount == 3 && (tree.getRoot()).isFake()) {
            return Double.NEGATIVE_INFINITY;
        }
        Node i;
        do {
			i = tree.getNode(Randomizer.nextInt(nodeCount));
        } while (i.isRoot() || i.getParent().isRoot() || (i).isDirectAncestor());

        final Node iParent = i.getParent();
        final Node iGrandParent = iParent.getParent();
        Node iUncle = iGrandParent.getLeft();
        if (iUncle.getNr() == iParent.getNr()) {
            iUncle = iGrandParent.getRight();
        }
        
        Boolean iIsRecipientBefore = Tools.isRecipient(i);
        Boolean iUncleIsRecipientBefore = Tools.isRecipient(iUncle);
        
        List<Integer> trTreeNodesNoTrOnGenes = Arrays.asList(Tools.getTrNodeNrsNotTransmissionOnGenes(tree,
				geneTreeIntervalsInput.get(), true));

        if (iUncle.getHeight() < iParent.getHeight()) {
			Tools.exchangeNodesRandomDirection(i, iUncle, iParent, iGrandParent);
			Boolean iIsRecipientAfter = Tools.isRecipient(i);
	        Boolean iUncleIsRecipientAfter = Tools.isRecipient(iUncle);
			if (trTreeNodesNoTrOnGenes.contains(iParent.getNr())) {
				if (!iIsRecipientBefore && iUncleIsRecipientAfter) {
					return Double.NEGATIVE_INFINITY;
				} else if (iIsRecipientBefore && !iUncleIsRecipientAfter) {
					return Double.NEGATIVE_INFINITY;
				}
			}
			if (trTreeNodesNoTrOnGenes.contains(iGrandParent.getNr())) {
				if (!iUncleIsRecipientBefore && iIsRecipientAfter) {
					return Double.NEGATIVE_INFINITY;
				} else if (iUncleIsRecipientBefore && !iIsRecipientAfter) {
					return Double.NEGATIVE_INFINITY;
				}
			}

			return 0.0;
        } else {
            // Couldn't find valid narrow move on this beast.tree!!
            return Double.NEGATIVE_INFINITY;
        }



    }

	/**
	 * @param tree
	 * @return log of hastings ratio for the wide move
	 */
	public double wide(final Tree tree) {

        final int nodeCount = tree.getNodeCount();

        //make sure that there are at least two distinct non-root nodes which are not direct ancestors.
        if (nodeCount == 3 && (tree.getRoot()).isFake()) {
            return Double.NEGATIVE_INFINITY;
        }

        Node i, j, iP, jP;
        do {
            i = tree.getNode(Randomizer.nextInt(nodeCount));
        } while (i.isRoot() || (i).isDirectAncestor());

        do {
            j = tree.getNode(Randomizer.nextInt(nodeCount));
        } while (j.getNr() == i.getNr() || j.isRoot() || (j).isDirectAncestor());

        iP = i.getParent();
        jP = j.getParent();

		Boolean iIsRecipientBefore = Tools.isRecipient(i);
		Boolean jIsRecipientBefore = Tools.isRecipient(j);

		Node CiP = getOtherChild(iP, i);
		Node CjP = getOtherChild(jP, j);

		List<Integer> trTreeNodesNoTrOnGenes = Arrays.asList(Tools.getTrNodeNrsNotTransmissionOnGenes(tree,
				geneTreeIntervalsInput.get(), false));

		if (trTreeNodesNoTrOnGenes.contains(i) || trTreeNodesNoTrOnGenes.contains(j)) {
			return Double.NEGATIVE_INFINITY;
		}

        if ((iP != jP) && (i != jP) && (j != iP)
                && (j.getHeight() < iP.getHeight())
                && (i.getHeight() < jP.getHeight())) {
			Tools.exchangeNodesRandomDirection(i, j, iP, jP);
			Boolean iIsRecipientAfter = Tools.isRecipient(i);
			Boolean jIsRecipientAfter = Tools.isRecipient(j);

			if (trTreeNodesNoTrOnGenes.contains(CiP.getNr())) {
				if (iIsRecipientBefore && !jIsRecipientAfter) {
					return Double.NEGATIVE_INFINITY;
				}
//				
//				if (!iIsRecipientBefore && jIsRecipientAfter) {
//					return Double.NEGATIVE_INFINITY;
//				} else if (iIsRecipientBefore && !jIsRecipientAfter) {
//					return Double.NEGATIVE_INFINITY;
//				}
			}
			if (trTreeNodesNoTrOnGenes.contains(CiP.getNr())) {
				if (jIsRecipientBefore && !iIsRecipientAfter) {
					return Double.NEGATIVE_INFINITY;
				}
//				if (!jIsRecipientBefore && iIsRecipientAfter) {
//					return Double.NEGATIVE_INFINITY;
//				} else if (jIsRecipientBefore && !iIsRecipientAfter) {
//					return Double.NEGATIVE_INFINITY;
//				}
			}

            return 0.0;
        }

        // Couldn't find valid wide move on this beast.tree!
        return Double.NEGATIVE_INFINITY;
    }


}
