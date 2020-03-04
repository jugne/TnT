package tnt.geneTree;

import java.io.PrintStream;

import beast.core.StateNode;
import beast.evolution.tree.TreeInterface;

public interface GeneTreeInterface extends TreeInterface {
	void init(PrintStream out);

	void close(PrintStream out);

	StateNode getCurrent();

}
