package tnt.transmissionTree;

import beast.evolution.tree.Node;
import beast.util.TreeParser;
import starbeast2.SpeciesTree;
import starbeast2.SpeciesTreeInterface;


public class TransmissionTree extends SpeciesTree implements SpeciesTreeInterface {


	/**
	 * reconstruct tree from XML fragment in the form of a DOM node *
	 */
	@Override
	public void fromXML(final org.w3c.dom.Node node) {
		final String newick = node.getTextContent();
		final TreeParser parser = new TreeParser();
		try {
			parser.thresholdInput.setValue(1e-10, parser);
		} catch (Exception e1) {
			e1.printStackTrace();
		}
		try {
			parser.offsetInput.setValue(0, parser);
			setRoot(parser.parseNewick(newick));
		} catch (Exception e) {
			e.printStackTrace();
		}

		initArrays();
		orientateTree();
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

	/**
	 * Orientates the tree according to stored (donor-recipient) metadata.
	 */
	public void orientateTree() {
		orientateNodeChildren(getRoot().getNr());
	}


	/**
	 * Adds the orientation (donor-recipient) metadata.
	 */
	public void addOrientationMetadata() {
		addOrientationMetadata(getRoot().getNr());
	}

	@Override
	public void store() {
		// add orientation metadata before storing tree
		// it is done since BEAST sort the tree in various places and changes the
		// orientation of child nodes
		addOrientationMetadata();
		super.store();
	}

	/**
	 * Orientate node children depending on stored metadata.
	 *
	 * @param subtreeRootNr the none number
	 */
	private void orientateNodeChildren(int subtreeRootNr) {
		Node subTreeRoot = this.getNode(subtreeRootNr);
		if (!subTreeRoot.isLeaf()) {
			if (!subTreeRoot.getLeft().metaDataString.equals("orientation=donor")) {
				Node left = subTreeRoot.getRight();
				Node right = subTreeRoot.getLeft();

				subTreeRoot.removeAllChildren(false);

				subTreeRoot.addChild(left);
				subTreeRoot.addChild(right);
			}

			orientateNodeChildren(subTreeRoot.getLeft().getNr());
			orientateNodeChildren(subTreeRoot.getRight().getNr());
		}

	}

	/**
	 * Add orientation metadata to each node. Left (0) child is always a donor.
	 * Right (1) child is always a recipient. Allows for: - tree state to be stored
	 * and restored from file even when BEAST applies sorting. - metadata can be
	 * used to color the output tree lineages for donors and recipients.
	 * 
	 * @param subtreeRootNr the none number
	 */
	private void addOrientationMetadata(int subtreeRootNr) {
		if (this.getNode(subtreeRootNr).isRoot()) {
			this.getNode(subtreeRootNr).metaDataString = "orientation=donor";
		}

		if (!this.getNode(subtreeRootNr).isLeaf()) {
			if (this.getNode(subtreeRootNr).isFake()) {
				this.getNode(subtreeRootNr).getLeft().metaDataString = this.getNode(subtreeRootNr).metaDataString;
				this.getNode(subtreeRootNr).getRight().metaDataString = this.getNode(subtreeRootNr).metaDataString;
			} else {
				this.getNode(subtreeRootNr).getLeft().metaDataString = "orientation=donor";
				this.getNode(subtreeRootNr).getRight().metaDataString = "orientation=recipient";
			}

			addOrientationMetadata(this.getNode(subtreeRootNr).getLeft().getNr());
			addOrientationMetadata(this.getNode(subtreeRootNr).getRight().getNr());
		}
	}
}
