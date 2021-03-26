package tnt.operators;

import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;
import starbeast2.SpeciesTreeInterface;
import tnt.distribution.GeneTreeIntervals;
import tnt.transmissionTree.TransmissionTree;
/**
 * @author Alexandra Gavryushkina
 */
@Description("")
public class TransmissionTreeScaleOperator extends beast.evolution.operators.ScaleOperator {

	public Input<SpeciesTreeInterface> transmissionTreeInput = new Input<>("transmissionTree",
			"Fully labeled transmission tree on which to simulate gene trees", Validate.REQUIRED);
	public Input<List<GeneTreeIntervals>> geneTreeIntervalsInput = new Input<>("geneTreeIntervals",
			"intervals for a gene tree", new ArrayList<>(), Validate.REQUIRED);

    //public Input<Boolean> m_pScaleSNodes = new Input<Boolean>("scaleSampledNodes", "If it is true then sampled node dates are scaled (default false).", false);

    @Override   //WARNING works with bifurcating (exactly 2 children) trees only
    // sampled ancestors are assumed to be on zero branches

    public double proposal() {

        final double scale = getScaler();
        final boolean scaleSNodes = false; // m_pScaleSNodes.get();

        try {

//            if () {
			TransmissionTree tree = (TransmissionTree) transmissionTreeInput.get();
                if (rootOnlyInput.get()) {
                    Node root = tree.getRoot();
                    if ((root).isFake() && !scaleSNodes) {
                        return Double.NEGATIVE_INFINITY;
                    }
                    double fNewHeight = root.getHeight() * scale;

                    //make sure the new height doesn't make a parent younger than a child
                    double oldestChildHeight;
                    if ((root).isFake()) {
                        oldestChildHeight = root.getNonDirectAncestorChild().getHeight();
                    } else oldestChildHeight = Math.max(root.getLeft().getHeight(), root.getRight().getHeight());
                    if (fNewHeight < oldestChildHeight) {
                        return Double.NEGATIVE_INFINITY;
                    }

                    root.setHeight(fNewHeight);

                    return -Math.log(scale);
                } else {
                    // scale the beast.tree
                    final int nScaledDimensions = tree.scale(scale);
                    //final int nScaledDimensions = tree.scale(scale, scaleSNodes);
                    return Math.log(scale) * (nScaledDimensions - 2);
                }
//            }
//            return Double.NEGATIVE_INFINITY;

        }  catch (Exception e) {
            return Double.NEGATIVE_INFINITY;
        }
    }

}
