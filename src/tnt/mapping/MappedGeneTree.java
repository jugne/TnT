package tnt.mapping;

import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.util.Log;
import beast.evolution.alignment.*;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import pitchfork.Pitchforks;
import tnt.distribution.GeneTreeIntervals;
import tnt.transmissionTree.TransmissionTree;
import tnt.util.MemoryFriendlyAlignment;
import tnt.util.MemoryFriendlyAncestralStateTreeLikelihood;
import tnt.util.Tools;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


/**
 * Maps hidden host switch events on a birth-death transmission tree.
 * 
 * @author Ugne Stolz <ugne.stolz@protonmail.com>
 * @date 7 Mar 2022
 */
public class MappedGeneTree extends Tree {

	public Input<TransmissionTree> trTreeInput = new Input<>("transmissionTree", "Transmission tree with NO hidden host switch events.",
			Input.Validate.REQUIRED);

	public Input<MappedTree> mappedTrTreeInput = new Input<>("mappedTransmissionTree", "Transmission tree WITH hidden host switch events.");

	public Input<Tree> geneTreeInput = new Input<>("geneTree", "Gene tree with NO hidden host switch events.");

	public Input<GeneTreeIntervals> geneTreeIntervalsInput = new Input<>("geneTreeIntervals",
			"intervals for a gene tree", Input.Validate.REQUIRED);

	final public Input<Alignment> dataInput = new Input<>("data", "sequence data for the beast.tree", Input.Validate.REQUIRED);

	public Input<SiteModel> siteModelInput = new Input<>("siteModel", "site model for sequence reconstruction");
	public Input<BranchRateModel> branchRateModelInput = new Input<>("branchRateModel", "branch rate model model for sequence reconstruction");

	public Input<String> outFileNameInput = new Input<>("outFileName",
			"Name of file to write sequences to.");

	public Input<Integer> burninInput = new Input<>("burnin",
			"How many states to skip before start logging.", 0);

	public Input<Boolean> logAllHiddenInput = new Input<>("logAllHidden",
			"Should all hidden sequences be logged or only ones informed by genetic data", false);

	TransmissionTree unmappedTrTree;
	MappedTree mappedTrTree;
	Tree unmappedGeneTree;
	Tree mappedGeneTree;
	GeneTreeIntervals intervals;
	Integer[] geneTreeNodeAssignment;

	static Integer[] trNodesToMappedNodes;

	int fakeLeafNr;
	Alignment data;
	int seqLength;
	List<Sequence> seqs = new ArrayList<>();
	MemoryFriendlyAlignment alignment = new MemoryFriendlyAlignment();

	MemoryFriendlyAncestralStateTreeLikelihood logger = new MemoryFriendlyAncestralStateTreeLikelihood();

	PrintStream ps;
	Boolean filtered = false;

	IntegerParameter constantsitesWeights;
	String filter;
	Boolean externalSeqFile = false;

	@Override
	public void initAndValidate() {
		super.initAndValidate();

		unmappedTrTree = trTreeInput.get();
		mappedTrTree = mappedTrTreeInput.get();
		unmappedGeneTree = geneTreeInput.get();
		intervals = geneTreeIntervalsInput.get();

		if (dataInput.get() instanceof FilteredAlignment){
			Log.warning.println("Only applicable when filter is for excluded constant sites! Do not use for any other filter yet!");
			data = ((FilteredAlignment) dataInput.get()).alignmentInput.get();
			constantsitesWeights = ((FilteredAlignment) dataInput.get()).constantSiteWeightsInput.get();
			filter = ((FilteredAlignment) dataInput.get()).filterInput.get();
			filtered = true;
		} else
			data = dataInput.get();



		assignFromWithoutID(new Tree(unmappedGeneTree.getRoot().copy()));
		if (outFileNameInput.get()!=null){
			externalSeqFile = true;
			try {
				ps = new PrintStream(outFileNameInput.get());
			} catch (FileNotFoundException e) {
				throw new RuntimeException("Error writing to output file '"
						+ outFileNameInput.get() + "'.");
			}
			ps.println("Sample"+"\t"+"Transmission Type"+"\t"+"Time"+"\t"+"Sequence");
		}

	}


	public void map() {

		unmappedTrTree = trTreeInput.get();
		if (mappedTrTreeInput.get() != null)
			mappedTrTree = mappedTrTreeInput.get();
		mappedGeneTree = geneTreeInput.get().copy();
		intervals = geneTreeIntervalsInput.get();

		if (mappedTrTreeInput.get()!=null){
			trNodesToMappedNodes = new Integer[unmappedTrTree.getNodeCount()];
			fillTrNodeMap(unmappedTrTree.getRoot(), mappedTrTree.getRoot());
		}

		geneTreeNodeAssignment = intervals.getGeneTreeNodeAssignment();
		fakeLeafNr = geneTreeInput.get().getLeafNodeCount();
		paintTransmissionNodes(mappedGeneTree.getRoot());


		// Ensure internal nodes are numbered correctly. (Leaf node numbers and
		// labels are matched to those in the untyped tree during the simulation.)
		numberInternalNodesOnSubtree(mappedGeneTree.getRoot(), fakeLeafNr);

		assignFromWithoutID(new Tree(mappedGeneTree.getRoot().copy()));

		seqs.clear();
		for (Sequence seq : data.sequenceInput.get()){
			Sequence nSeq = new Sequence(seq.getTaxon(), seq.getData());
			seqs.add(nSeq);
		}

		seqLength = seqs.get(0).dataInput.get().length();
		char[] charArray = new char[seqLength];
		Arrays.fill(charArray, '-');
		String newString = new String(charArray);
		for (int i=geneTreeInput.get().getLeafNodeCount(); i<fakeLeafNr; i++)
			seqs.add(new Sequence(Integer.toString(i), newString));

		if (filtered){
			Alignment al = new Alignment();
			al.initByName("sequence", seqs, "statecount", 4);
			alignment.initByName("constantSiteWeights", constantsitesWeights, "filter", filter,
					"data", al);
		} else
			alignment.initByName("sequence", seqs, "statecount", 4);

//if (first) {
//		logger = new AncestralStateTreeLikelihood();
		logger.initByName("data", alignment, "siteModel", siteModelInput.get(),"branchRateModel", branchRateModelInput.get(), "tree", this, "tag", "seq", "sampleTips", false);
//first = false;
//}

		logger.calculateLogP();
		logger.redrawAncestralStates();

	}

	// geneRoot has to be a logicalNode
	private void paintTransmissionNodes(Node geneRoot){
		Node geneNode = geneRoot;

		int trTreeNr = geneTreeNodeAssignment[geneNode.getNr()];
		Node TrTreeNode = unmappedTrTree.getNode(trTreeNr);
		if (mappedTrTreeInput.get()!=null)
			TrTreeNode = mappedTrTree.getNode(trNodesToMappedNodes[trTreeNr]);
		boolean recipient = TrTreeNode.metaDataString.contains("recipient");

		if (!geneNode.isRoot()){
			Node logicalParent = Pitchforks.getLogicalParent(geneNode);
			if (TrTreeNode.isLeaf())
				TrTreeNode = TrTreeNode.getParent();
			while (Tools.greaterHeightNode(logicalParent, TrTreeNode)){
				if (!TrTreeNode.isFake() && (TrTreeNode.getChildCount()==1 || recipient) &&
						Tools.greaterHeightNode(TrTreeNode, geneNode)){
					Node tmp = new Node();
					Node parent = geneNode.getParent();
					parent.removeChild(geneNode);
					parent.addChild(tmp);
					tmp.addChild(geneNode);
					Node tmpChild = new Node();
					tmpChild.setNr(fakeLeafNr);
					tmpChild.setID(Integer.toString(fakeLeafNr));
					fakeLeafNr +=1;
					tmp.addChild(tmpChild);
					tmp.setHeight(TrTreeNode.getHeight());
					tmpChild.setHeight(TrTreeNode.getHeight());
					if (TrTreeNode.getChildCount()==1)
						tmp.metaDataString="nodeType=hidden";
					else if (recipient)
						tmp.metaDataString="nodeType=observed_l";
					geneNode = tmp;
				}
				recipient = TrTreeNode.metaDataString.contains("recipient");
				if (TrTreeNode.isRoot())
					break;
				TrTreeNode = TrTreeNode.getParent();
			}
			if (logicalParent.metaDataString==null && !TrTreeNode.isLeaf()
					&& !TrTreeNode.isFake())
				if(Tools.equalHeightWithPrecision(TrTreeNode, logicalParent)) {
					logicalParent.metaDataString = "nodeType=observed";
				} else if(Tools.isMultiMerger(Tools.getLogicalNotHiddenChildren(mappedGeneTree.getRoot()),
						logicalParent) || Pitchforks.isPolytomy(logicalParent)){
					logicalParent.metaDataString = "nodeType=hidden_gene";
				}
			}

		List<Node> children = Pitchforks.getLogicalChildren(geneNode);
		for (Node child : children){
			paintTransmissionNodes(child);
		}


	}

	private static void fillTrNodeMap(Node transmissionSubRoot, Node mappedTrSubRoot){
		while (mappedTrSubRoot.getChildren().size()==1){
			mappedTrSubRoot = mappedTrSubRoot.getChild(0);
		}
		trNodesToMappedNodes[transmissionSubRoot.getNr()] = mappedTrSubRoot.getNr();
		if (!transmissionSubRoot.isLeaf()) {
			fillTrNodeMap(transmissionSubRoot.getChild(0), mappedTrSubRoot.getChild(0));
			fillTrNodeMap(transmissionSubRoot.getChild(1), mappedTrSubRoot.getChild(1));
		}
	}

	/**
	 * COPIED from BDMM prime package by T. G. Vaughan on 2022-03-12
	 *
	 * Apply node numbers to internal nodes below and including subtreeRoot. Numbers
	 * are applied postorder, so parents always have larger numbers than their
	 * children and the root has the hightest number.
	 *
	 * @param subtreeRoot root of subtree
	 * @param nextNumber  next number to be used
	 * @return next number to be used on another part of the tree.
	 */
	private int numberInternalNodesOnSubtree(Node subtreeRoot, int nextNumber) {

		if (subtreeRoot.isLeaf())
			return nextNumber;

		for (Node child : subtreeRoot.getChildren())
			nextNumber = numberInternalNodesOnSubtree(child, nextNumber);

		subtreeRoot.setNr(nextNumber);

		return nextNumber + 1;
	}

	private long lastRemapSample = -1;

	/**
	 * Remap the tree.  Intended to be called by loggers requiring
	 * a mapped tree.  Supplying the sample number allows the result
	 * of the remapping to be cached and used for other loggers.
	 *
	 * @param sample sample number at log
	 */
	public void remapForLog(long sample) {
		if (sample == lastRemapSample)
			return;

		map();
		lastRemapSample = sample;
	}

	@Override
	public void init(PrintStream out) {
		if (!externalSeqFile)
			unmappedGeneTree.init(out);
	}

	@Override
	public void log(long sample, PrintStream out) {
		if (sample > burninInput.get()) {
			remapForLog(sample);

			Tree tree = (Tree) getCurrent();
			if (externalSeqFile) {
				logSeqsToFile(sample, tree.getRoot());
				seqs.clear();
				alignment.finalize();
				try {
					logger.finalize();
				} catch (Throwable e) {
					throw new RuntimeException(e);
				}
			} else {
				// write out the log tree with metadata
				out.print("tree STATE_" + sample + " = ");
				out.print(toNewick(tree.getRoot()));
				out.print(";");
			}
		}
	}

	@Override
	public void close(PrintStream out) {
		if (externalSeqFile)
			geneTreeInput.get().close(out);
	}

	void logSeqsToFile(long sample, Node node){
		if (node.metaDataString != null){
			if (logAllHiddenInput.get() ||
					(node.metaDataString.contains("gene") ||
							node.metaDataString.contains("observed"))) {
				int[] patternStates = logger.getStatesForNode(this, node);
				int[] siteStates = new int[seqLength];
				for (int i = 0; i < seqLength; i++) {
					siteStates[i] = patternStates[alignment.getPatternIndex(i)];
				}
				StringBuffer buf = new StringBuffer();
				String seq = logger.getDataType().encodingToString(siteStates);
				buf.append(sample);
				buf.append("\t");
				buf.append(node.metaDataString.split("=")[1]);
				buf.append("\t");
				buf.append(node.getHeight());
				buf.append("\t");
				buf.append(node.getNr());
				buf.append("\t");
				buf.append(seq);

				ps.println(buf);
				ps.flush();
			}
		}
		for (Node child : node.getChildren()){
			logSeqsToFile(sample, child);
		}
	}

	String toNewick(Node node) {
		StringBuilder buf = new StringBuilder();
		if (node.getLeft() != null) {
			buf.append("(");
			buf.append(toNewick(node.getLeft()));
			if (node.getRight() != null) {
				buf.append(',');
				buf.append(toNewick(node.getRight()));
			}
			buf.append(")");
		}
		buf.append(node.getNr() + 1);

		int [] patternstates = logger.getStatesForNode(this, node);
		int[] siteStates = new int[seqLength];
		for (int i = 0; i < seqLength; i++) {
			siteStates[i] = patternstates[alignment.getPatternIndex(i)];
		}
		String seq = logger.getDataType().encodingToString(siteStates);

		if (node.metaDataString != null){
				buf.append("[&");
				buf.append(node.metaDataString);
				buf.append(",");
				buf.append("nr" + "=").append(node.getNr());
				buf.append(",");
				buf.append("seq" + "=\"").append(seq).append("\"");
				buf.append(']');
		}


		buf.append(":");

		double nodeLength;
		nodeLength = node.getLength();

		buf.append(nodeLength);

		return buf.toString();
	}

}
