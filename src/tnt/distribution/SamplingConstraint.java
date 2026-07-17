package tnt.distribution;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import tnt.mapping.MappedTree;
import tnt.util.Tools;

public class SamplingConstraint extends Distribution {

	public Input<Tree> trTreeInput = new Input<>("transmissionTree",
			"transmission tree");
	public Input<Tree> treeInput = new Input<>("tree", "tree input", Input.Validate.XOR, trTreeInput);
	public Input<List<TaxonSet>> taxonsetsInput = new Input<>("taxonsets",
			"a separate list of taxa for samples collected from the same patient", new ArrayList<>());

	public Input<MappedTree> mapTrTreeInput = new Input<>("mappedTransmissionTree",
			"mapped transmission tree");
	public Input<TaxonSet> mrcaTaxonSetInput = new Input<>("mrcaTaxonSet",
			"a separate list of taxa which should be direct descendants of the root, " +
					"should be only used with mapped tree");

	public Input<Integer> maxMappedNodesInput = new Input<>("maxMappedNodes",
			"maximum unobserved transmission events, " +
					"should be only used with mapped tree", Integer.MAX_VALUE);



	private List<List<Node>> taxonsets;
	private List<Node> mrcaTaxonSet;
	private List<String> mrcaTaxonSetHostList;
	private int nSets;
	private int nSamples;

	private boolean isStochastic = false;



	@Override
	public void initAndValidate() {
		Tree trTree = trTreeInput.get();
		if (trTree == null)
			trTree = treeInput.get();
		MappedTree mapTrTree = mapTrTreeInput.get();
		if (mapTrTree != null){
			mapTrTree.addOrientationMetadata();
			mapTrTree.addHostMetadata();
		}


		nSets = taxonsetsInput.get().size();
		nSamples = trTree.getLeafNodeCount();
		taxonsets = new ArrayList<List<Node>>(nSets);
		if (mrcaTaxonSetInput.get() !=null)
			mrcaTaxonSet = new ArrayList<Node>(mrcaTaxonSetInput.get().getTaxonCount());

		for (int i = 0; i < nSets; i++)
			taxonsets.add(new ArrayList<Node>());


		for (int idx = 0; idx < nSamples; idx++) {
			Node leaf = trTree.getNode(idx);
			int i = 0;
			for (TaxonSet t : taxonsetsInput.get()) {
				List<String> tmp = t.asStringList();
				if (tmp.contains(leaf.getID())) {
					taxonsets.get(i).add(leaf);
					break;
				} else {
					i += 1;
				}
			}

		}
		for (int ii = 0; ii < nSets; ii++) {
			taxonsets.get(ii).sort(Comparator.comparing(Node::getHeight));
		}

		if(mapTrTree != null && mrcaTaxonSetInput.get() != null){
			for (int idx=0; idx < nSamples; idx++){
				Node leaf = mapTrTree.getNode(idx);
				List<String> tmp = mrcaTaxonSetInput.get().asStringList();
				if (tmp.contains(leaf.getID())) {
					mrcaTaxonSet.add(leaf);
				}
			}



			mrcaTaxonSetHostList = new ArrayList<String>();

			for (int ii = 0; ii < mrcaTaxonSet.size(); ii++) {
				mrcaTaxonSetHostList.add(mrcaTaxonSet.get(ii).getID());
			}
			mrcaTaxonSetHostList.add("unsampled");
			mrcaTaxonSet.sort(Comparator.comparing(Node::getHeight));
			if (mrcaTaxonSetHostList.size()>0)
				isStochastic = true;
		}
	}

	@Override
	public boolean isStochastic() {
		return isStochastic;
	}

	@Override
	public double calculateLogP() {
		logP = 0;
		Tree trTree = trTreeInput.get();
		if (trTree == null)
			trTree = treeInput.get();
		for (int j = 0; j < nSets; j++) {
			List<Node> tmp = new ArrayList<>(taxonsets.get(j));
			Node startLeaf = tmp.get(0);
			Node child = trTree.getNode(startLeaf.getNr());
			Node parent = child.getParent();
			if (parent.getLeft().getNr() == child.getNr() && !parent.getLeft().isDirectAncestor()) {
				tmp.remove(startLeaf);
			}
			while (parent != null) {
				if (parent.isFake()) {
					if (tmp.get(0).getNr() == parent.getDirectAncestorChild().getNr()) {
						if (!tmp.remove(tmp.get(0)))
							System.out.println("SamplingThroughTime distribution error");
						if (tmp.isEmpty())
							break;
					} else {
						logP = Double.NEGATIVE_INFINITY;
						return logP;
					}
				} else if (!parent.isRoot() && parent.getChild(0).getNr() != child.getNr()) {
					logP = Double.NEGATIVE_INFINITY;
					return logP;
				}
				child = parent;
				parent = child.getParent();
			}
			if (!tmp.isEmpty()) {
				logP = Double.NEGATIVE_INFINITY;
				return logP;
			}

		}

		if (mrcaTaxonSet!=null){
			int i=0;
			if (mrcaTaxonSet.size() != 0) {
				if (!mapTrTreeInput.get().map()){
					logP = Double.NEGATIVE_INFINITY;
					return logP;
				}
				if(mapTrTreeInput.get().getTotalNodesNumber() > maxMappedNodesInput.get()){
					logP = Double.NEGATIVE_INFINITY;
					return logP;
				}
				mapTrTreeInput.get().addOrientationMetadata();
				mapTrTreeInput.get().addHostMetadata();
				mapTrTreeInput.get().setEverythingDirty(true);

//			Node nn = mapTrTreeInput.get().getRoot();
//			while(!nn.isLeaf()){
//				if (nn.getLeft().getMetaData("host").equals("unsampled"))
//					System.out.println();
//				nn = nn.getLeft();
//			}

			}

//		System.out.println(mapTrTreeInput.get().getRoot().toNewick());

			for (Node l : mrcaTaxonSet){
				Node parent = mapTrTreeInput.get().getNode(l.getNr()).getParent();
				while (parent != null) {
//				if (parent.isRoot()) {
//					break;
//				}
//				if (parent.getMetaData("host").toString().equals("unsampled")) {
//					System.out.println();
//				}
//				if (!mrcaTaxonSetHostList.contains(parent.getMetaData("host").toString())) {
					if (!mrcaTaxonSetHostList.contains(Tools.extractHostValue(parent.metaDataString))) {
						i++;
						break;
					}
					parent = parent.getParent();
				}
			}
			if(i==mrcaTaxonSet.size()){
				logP = Double.NEGATIVE_INFINITY;
				return logP;
			}
		}



		return 0.0;
	}


	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub

	}

}
