package tnt.operators;

import java.util.ArrayList;
import java.util.List;

import beast.core.Input;
import beast.core.Input.Validate;
import starbeast2.SpeciesTreeInterface;
import tnt.distribution.GeneTreeIntervals;

public class AllTreeScale extends beast.evolution.operators.ScaleOperator {
	public Input<SpeciesTreeInterface> transmissionTreeInput = new Input<>("transmissionTree",
			"Fully labeled transmission tree on which to simulate gene trees", Validate.REQUIRED);
	public Input<List<GeneTreeIntervals>> geneTreeIntervalsInput = new Input<>("geneTreeIntervals",
			"intervals for a gene tree", new ArrayList<>(), Validate.REQUIRED);

	pitchfork.operators.ScaleOperator treeScale;

	@Override
	public void initAndValidate() {
		super.initAndValidate();
	}

}
