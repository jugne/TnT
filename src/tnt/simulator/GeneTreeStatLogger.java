package tnt.simulator;



import java.io.PrintStream;

import beast.core.Description;
import beast.core.Function;
import beast.core.Loggable;
import beast.evolution.tree.TreeStatLogger;


@Description("Logger to report statistics of a tree")
public class GeneTreeStatLogger extends TreeStatLogger implements Loggable, Function {

    @Override
    public void initAndValidate() {
		super.initAndValidate();

    }

    @Override
    public void init(PrintStream out) {
		super.init(out);
		out.print("hiddenNodeCount\t");
		out.print("hiddenNodeTimes\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
		super.log(sample, out);
		final SimulatedGeneTree tree = (SimulatedGeneTree) super.treeInput.get();
		out.print(tree.hiddenNodes + "\t");
		out.print(tree.hiddenNodeTimes + "\t");
    }
}
