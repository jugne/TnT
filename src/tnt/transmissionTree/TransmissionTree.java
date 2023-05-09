package tnt.transmissionTree;

import beast.evolution.tree.Node;
import beast.util.TreeParser;
import starbeast2.SpeciesTree;
import starbeast2.SpeciesTreeInterface;
import tnt.util.Tools;


public class TransmissionTree extends SpeciesTree implements SpeciesTreeInterface {

	private boolean hasOrientationMetadata = false;
	private boolean hasHostMetadata = false;


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

	public boolean hasOrientationMetadata(){
		if (this.getRoot().metaDataString!=null && this.getRoot().metaDataString.contains("orientation")){
			hasOrientationMetadata = true;
			return true;
		}
		hasOrientationMetadata = false;
		return false;
	}
	public boolean hasHostMetadata(){
		if (this.getRoot().metaDataString!=null && this.getRoot().metaDataString.contains("host")){
			hasHostMetadata = true;
			return true;
		}
		hasHostMetadata = false;
		return false;
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
		hasOrientationMetadata = true;
	}

	public void addHostMetadata() {
		if (!hasOrientationMetadata()){
			addOrientationMetadata();
		}
		addHostMetadata(getRoot().getNr());
		hasHostMetadata = true;
	}

	@Override
	public void store() {
		// add orientation metadata before storing tree
		// it is done since BEAST sort the tree in various places and changes the
		// orientation of child nodes
//		if (!hasHostMetadata) {
			addOrientationMetadata();
//		}
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
			if ((!subTreeRoot.isFake() && !subTreeRoot.getChild(0).metaDataString.contains("orientation=donor"))
					|| (subTreeRoot.isFake() && subTreeRoot.getChild(1).getHeight() != subTreeRoot.getHeight())) {
				Node left = subTreeRoot.getChild(1);
				Node right = subTreeRoot.getChild(0);

				subTreeRoot.removeAllChildren(false);

				subTreeRoot.setLeft(left);
				subTreeRoot.setRight(right);
//				subTreeRoot.addChild(left);
//				subTreeRoot.addChild(right);
			}

			orientateNodeChildren(subTreeRoot.getChild(0).getNr());
			orientateNodeChildren(subTreeRoot.getChild(1).getNr());
		}
	}

	/**
	 * Add orientation metadata to each node. Left (0) child is always a donor.
	 * Right (1) child is always a recipient. Allows for: - tree state to be stored
	 * and restored from file even when BEAST applies sorting. - metadata can be
	 * used to color the output tree lineages for donors and recipients.
	 * 
	 * @param subtreeRootNr the node number
	 */
	private void addOrientationMetadata(int subtreeRootNr) {
		Node subRoot = this.getNode(subtreeRootNr);
		if (subRoot.isRoot()) {
			subRoot.metaDataString = "orientation=donor";
		}

		if (!subRoot.isLeaf()) {
			if (subRoot.isFake()) {
				subRoot.getLeft().metaDataString = subRoot.metaDataString;
				subRoot.getRight().metaDataString = subRoot.metaDataString;
			} else if (subRoot.getChildCount()==1){
				subRoot.getLeft().metaDataString = subRoot.metaDataString;
			} else {
				subRoot.getLeft().metaDataString = "orientation=donor";
				subRoot.getRight().metaDataString = "orientation=recipient";
			}

			addOrientationMetadata(subRoot.getLeft().getNr());
			if(subRoot.getChildCount()!=1){
				addOrientationMetadata(subRoot.getRight().getNr());
			}
		}
	}

	/**
	 * Adds metadata on host occupying transmission tree edge.
	 * Should be used on trees with full transmission histories.
	 * That means either fully sampled trees or trees with stochastically mapped hidden transmission events.
	 *
	 * Relies on orientation metadata being added before
	 * (public parent mathod makes sure of this).
	 *
	 * @param subtreeRootNr the node number
	 */
	private void addHostMetadata(int subtreeRootNr){
		Node subTreeRoot = this.getNode(subtreeRootNr);
		String metaData = subTreeRoot.metaDataString;

		if (subTreeRoot.isLeaf()) {
			subTreeRoot.metaDataString = String.format("%s,%s=%s",
					metaData, "host", subTreeRoot.getID().split("_")[0]);
			return;
		}
		if (subTreeRoot.getChildCount()==1 && !subTreeRoot.isRoot()){
			Node parent = subTreeRoot.getParent();
			if (parent.getChildCount()>1 && !parent.isFake()){ // coalescent (transmission event)
				if (subTreeRoot.metaDataString.contains("orientation=donor")){
					subTreeRoot.metaDataString = parent.metaDataString;
				} else {
					subTreeRoot.metaDataString = String.format("%s,%s=%s",
							metaData, "host", "unsampled");
				}
			} else if(parent.isFake()){ // fake (sampling event)
				subTreeRoot.metaDataString = String.format("%s,%s=%s",
						metaData, "host", "unsampled");//parent.metaDataString;
			}
			else if (parent.getChildCount()==1){ // jump (unobserved transmission event)
				subTreeRoot.metaDataString = String.format("%s,%s=%s",
						metaData, "host", "unsampled");
			}
		} else if (subTreeRoot.getChildCount()==1 && subTreeRoot.isRoot()){
			subTreeRoot.metaDataString = String.format("%s,%s=%s",
					metaData, "host", "unsampled");
		} else if (subTreeRoot.isFake()){
			subTreeRoot.metaDataString = String.format("%s,%s=%s",
					metaData, "host", firstLeafId(subTreeRoot.getDirectAncestorChild()));
		}else {
			subTreeRoot.metaDataString = String.format("%s,%s=%s",
					metaData, "host", firstLeafId(subTreeRoot.getLeft())); // bug here
		}
			addHostMetadata(subTreeRoot.getLeft().getNr());
			if(subTreeRoot.getChildCount()!=1){
				addHostMetadata(subTreeRoot.getRight().getNr());
			}
	}

	private String firstLeafId(Node n) {

		if (n.isLeaf()) {
			return n.getID().split("_")[0];
		} else if (Tools.equalHeightWithPrecision(n, n.getLeft())) {
			return n.getLeft().getID().split("_")[0];
		} else if (n.getChildCount()>1 && Tools.equalHeightWithPrecision(n, n.getRight())) {
			return n.getRight().getID().split("_")[0];
		}else if (n.getChildCount()==1){
			return "unsampled";
		} else {
			return n.getChild(0).metaDataString.contains("orientation=donor") ? firstLeafId(n.getChild(0))
					: firstLeafId(n.getChild(1));
		}

	}
}
