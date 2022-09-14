package tnt.operators;

import java.util.ArrayList;
import java.util.List;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import tnt.distribution.GeneTreeIntervals;
import tnt.util.Tools;

public class SwapOrientations extends TreeOperator {

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

		// get fit nodes as list
		List<Integer> fitNodesNrs = new ArrayList<Integer>() {
			{
				for (int i : Tools.getTrNodeNrsNotTransmissionOnGenes(tree, geneTreeIntervalsInput.get(), true))
					add(i);
			}
		};

		List<Integer> tmp = new ArrayList<Integer>();
		for (int i : fitNodesNrs) {
			if (tree.getNode(i).isRoot() || tree.getNode(i).getParent().isFake())
				tmp.add(i);
		}
		fitNodesNrs.removeAll(tmp);


		if (fitNodesNrs.size() == 0)
			return Double.NEGATIVE_INFINITY;
		if (fitNodesNrs.size() < 0) {
			System.err.print("Wrong node counting in Swap operator!");
			System.exit(1);
		}

        Node i;

//        do {
		int nr = fitNodesNrs.get(Randomizer.nextInt(fitNodesNrs.size()));
		i = tree.getNode(nr);
//		} while (i.isRoot() || i.isDirectAncestor() || i.getParent().isFake());


        Node iP = i.getParent();
		Node left = iP.getLeft();
		Node right = iP.getRight();

		iP.setLeft(right);
		iP.setRight(left);
/*        Node CiP;
        if (iP.getLeft().getNr() == i.getNr()) {
            CiP = iP.getRight();
			iP.removeAllChildren(true);
			iP.addChild(CiP);
			iP.addChild(i);
        } else {
            CiP = iP.getLeft();
			iP.removeAllChildren(true);
			iP.addChild(i);
			iP.addChild(CiP);
        }*/


		return 0.0;

    }


}
