package tnt.transmissionTree;

import starbeast2.SpeciesTree;
import starbeast2.SpeciesTreeInterface;


public class TransmissionTree extends SpeciesTree implements SpeciesTreeInterface {


	/**
	 * reconstruct tree from XML fragment in the form of a DOM node *
	 */
	@Override
	public void fromXML(final org.w3c.dom.Node node) {
		final String newick = node.getTextContent();
		final TransmissionTreeParser parser = new TransmissionTreeParser();
		try {
			parser.thresholdInput.setValue(1e-10, parser);
		} catch (Exception e1) {
			e1.printStackTrace();
		}
		try {
			parser.offsetInput.setValue(0, parser);
			setRoot(parser.parseNewick(newick));
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		initArrays();
	}

	@Override
	public TransmissionTree copy() {
		TransmissionTree tree = new TransmissionTree();
		tree.setID(getID());
		tree.index = index;
		tree.root = root.copy();
		tree.nodeCount = nodeCount;
		tree.internalNodeCount = internalNodeCount;
		tree.leafNodeCount = leafNodeCount;
		return tree;
	}
}
