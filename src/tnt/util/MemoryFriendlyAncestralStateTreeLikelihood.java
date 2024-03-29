package tnt.util;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.UserDataType;
import beast.evolution.likelihood.LeafTrait;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.*;
import beast.util.Randomizer;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Alexei Drummond
 * @author Ugne Stolz
 */
@Description("Ancestral State Tree Likelihood. " +
        "Copy of AncestralStateTreeLikelihood from beast-classic package. " +
        "Only finalize method added to solve java memory leaks." +
        "No beagle support.")
public class MemoryFriendlyAncestralStateTreeLikelihood extends TreeLikelihood implements TreeTraitProvider {
    public static final String STATES_KEY = "states";

    public Input<String> tagInput = new Input<>("tag","label used to report trait", Validate.REQUIRED);
    public Input<Boolean> useMAPInput = new Input<>("useMAP","whether to use maximum aposteriori assignments or sample", false);
    public Input<Boolean> returnMLInput = new Input<>("returnML", "report integrate likelihood of tip data", true);

    public Input<Boolean> useJava = new Input<>("useJava", "prefer java, even if beagle is available", true);

    public Input<Boolean> sampleTipsInput = new Input<>("sampleTips", "if tips have missing data/ambigous values sample them for logging (default true)", true);

	public Input<List<LeafTrait>> leafTriatsInput = new Input<>("leaftrait", "list of leaf traits",
			new ArrayList<>());

    protected DataType dataType;
    private int[][] reconstructedStates;
    private int[][] storedReconstructedStates;
    private String tag;
    private boolean areStatesRedrawn = false;
    private boolean storedAreStatesRedrawn = false;
    private boolean useMAP = false;
    private boolean returnMarginalLogLikelihood = true;
    private double jointLogLikelihood;
    private double storedJointLogLikelihood;
    boolean likelihoodKnown = false;

	int[][] storedTipStates;

	/** parameters for each of the leafs **/
	IntegerParameter[] parameters;

	/** and node number associated with parameter **/
	int[] leafNr;
	int traitDimension;
    int patternCount;
    int stateCount;

    int[][] tipStates; // used to store tip states when using beagle

    @Override
    public void initAndValidate() {
    	if (dataInput.get().getSiteCount() == 0) {
    		return;
    	}


    	String sJavaOnly = null;
    	if (useJava.get()) {
    		sJavaOnly = System.getProperty("java.only");
    		System.setProperty("java.only", "" + true);
    	}
    	super.initAndValidate();
    	if (useJava.get()) {
	    	if (sJavaOnly != null) {
	    		System.setProperty("java.only", sJavaOnly);
	    	} else {
	    		System.clearProperty("java.only");
	    	}
    	}

        this.tag = tagInput.get();
        TreeInterface treeModel = treeInput.get();
        patternCount = dataInput.get().getPatternCount();
        dataType = dataInput.get().getDataType();
        stateCount = dataType.getStateCount();

        reconstructedStates = new int[treeModel.getNodeCount()][patternCount];
        storedReconstructedStates = new int[treeModel.getNodeCount()][patternCount];

        this.useMAP = useMAPInput.get();
        this.returnMarginalLogLikelihood = returnMLInput.get();

        treeTraits.addTrait(STATES_KEY, new TreeTrait.IA() {
            public String getTraitName() {
                return tag;
            }

            public Intent getIntent() {
                return Intent.NODE;
            }

            public int[] getTrait(TreeInterface tree, Node node) {
                return getStatesForNode(tree,node);
            }

            public String getTraitString(TreeInterface tree, Node node) {
                return formattedState(getStatesForNode(tree,node), dataType);
            }
        });

        if (beagle != null) {
            if (!(siteModelInput.get() instanceof SiteModel.Base)) {
            	throw new IllegalArgumentException ("siteModel input should be of type SiteModel.Base");
            }
            m_siteModel = (SiteModel.Base) siteModelInput.get();
        	substitutionModel = m_siteModel.substModelInput.get();
            int nStateCount = dataInput.get().getMaxStateCount();
            probabilities = new double[(nStateCount + 1) * (nStateCount + 1)];
        }

        int tipCount = treeModel.getLeafNodeCount();
        tipStates = new int[tipCount][];

        Alignment data = dataInput.get();
        for (Node node : treeInput.get().getExternalNodes()) {
            String taxon = node.getID();
            int taxonIndex = data.getTaxonIndex(taxon);
            if (taxonIndex == -1) {
            	if (taxon.startsWith("'") || taxon.startsWith("\"")) {
                    taxonIndex = data.getTaxonIndex(taxon.substring(1, taxon.length() - 1));
                }
                if (taxonIndex == -1) {
                	throw new RuntimeException("Could not find sequence " + taxon + " in the alignment");
                }
            }
            tipStates[node.getNr()] = new int[patternCount];
            if (!m_useAmbiguities.get()) {
            	likelihoodCore.getNodeStates(node.getNr(), tipStates[node.getNr()]);
            } else {
            	int [] states = tipStates[node.getNr()];
	            for (int i = 0; i < patternCount; i++) {
	                int code = data.getPattern(taxonIndex, i);
	                int[] statesForCode = data.getDataType().getStatesForCode(code);
	                if (statesForCode.length==1)
	                    states[i] = statesForCode[0];
	                else
	                    states[i] = code; // Causes ambiguous states to be ignored.
	            }
            }
    	}

        if (m_siteModel.getCategoryCount() > 1)
            throw new RuntimeException("Reconstruction not implemented for multiple categories yet.");


        // stuff for dealing with ambiguities in tips
        if (!m_useAmbiguities.get() && leafTriatsInput.get().size() == 0) {
        	return;
        }
		traitDimension = tipStates[0].length;

		leafNr = new int[leafTriatsInput.get().size()];
		parameters = new IntegerParameter[leafTriatsInput.get().size()];

		List<String> taxaNames = dataInput.get().getTaxaNames();
		for (int i = 0; i < leafNr.length; i++) {
			LeafTrait leafTrait = leafTriatsInput.get().get(i);
			parameters[i] = leafTrait.parameter.get();
			// sanity check
			if (parameters[i].getDimension() != traitDimension) {
				throw new IllegalArgumentException("Expected parameter dimension to be " + traitDimension + ", not "
						+ parameters[i].getDimension());
			}
			// identify node
			String taxon = leafTrait.taxonName.get();
			int k = 0;
			while (k < taxaNames.size() && !taxaNames.get(k).equals(taxon)) {
				k++;
			}
			leafNr[i] = k;
			// sanity check
			if (k == taxaNames.size()) {
				throw new IllegalArgumentException("Could not find taxon '" + taxon + "' in tree");
			}
			// initialise parameter value from states
			Integer[] values = new Integer[tipStates[k].length];
			for (int j = 0; j < tipStates[k].length; j++) {
				values[j] = tipStates[k][j];
			}
			IntegerParameter p = new IntegerParameter(values);
			p.setLower(0);
			p.setUpper(dataType.getStateCount()-1);
			parameters[i].assignFromWithoutID(p);
		}

		storedTipStates = new int[tipStates.length][traitDimension];
		for (int i = 0; i < tipStates.length; i++) {
			System.arraycopy(tipStates[i], 0, storedTipStates[i], 0, traitDimension);
		}
    }

    @Override
    public void store() {
        super.store();

        for (int i = 0; i < reconstructedStates.length; i++) {
            System.arraycopy(reconstructedStates[i], 0, storedReconstructedStates[i], 0, reconstructedStates[i].length);
        }

        storedAreStatesRedrawn = areStatesRedrawn;
        storedJointLogLikelihood = jointLogLikelihood;


        // deal with ambiguous tips
        if (leafNr != null) {
            for (int k : leafNr) {
                System.arraycopy(tipStates[k], 0, storedTipStates[k], 0, traitDimension);
            }
        }
    }

    @Override
    public void restore() {

        super.restore();

        int[][] temp = reconstructedStates;
        reconstructedStates = storedReconstructedStates;
        storedReconstructedStates = temp;

        areStatesRedrawn = storedAreStatesRedrawn;
        jointLogLikelihood = storedJointLogLikelihood;

        // deal with ambiguous tips
        if (leafNr != null) {
            for (int k : leafNr) {
                int[] tmp = tipStates[k];
                tipStates[k] = storedTipStates[k];
                storedTipStates[k] = tmp;
                // Does not handle ambiguities or missing taxa
                likelihoodCore.setNodeStates(k, tipStates[k]);
            }
        }
    }

    @Override
    protected boolean requiresRecalculation() {
    	likelihoodKnown = false;

    	boolean isDirty = super.requiresRecalculation();
    	if (!m_useAmbiguities.get()) {
    		return isDirty;
    	}

    	int hasDirt = Tree.IS_CLEAN;

		// check whether any of the leaf trait parameters changed
		for (int i = 0; i < leafNr.length; i++) {
			if (parameters[i].somethingIsDirty()) {
				int k = leafNr[i];
				for (int j = 0; j < traitDimension; j++) {
					tipStates[k][j] = parameters[i].getValue(j);
				}
				likelihoodCore.setNodeStates(k, tipStates[k]);
				isDirty = true;
				// mark leaf's parent node as dirty
				Node leaf = treeInput.get().getNode(k);
				// leaf.makeDirty(Tree.IS_DIRTY);
				leaf.getParent().makeDirty(Tree.IS_DIRTY);
	            hasDirt = Tree.IS_DIRTY;
			}
		}
		isDirty |= super.requiresRecalculation();
		this.hasDirt |= hasDirt;

		return isDirty;


    }


    public DataType getDataType() {
        return dataType;
    }

    public int[] getStatesForNode(TreeInterface tree, Node node) {
        if (tree != treeInput.get()) {
            throw new RuntimeException("Can only reconstruct states on treeModel given to constructor");
        }

        if (!likelihoodKnown) {
        	try {
        		 calculateLogP();
        	} catch (Exception e) {
				throw new RuntimeException(e.getMessage());
			}
        }

        if (!areStatesRedrawn) {
            redrawAncestralStates();
        }
        return reconstructedStates[node.getNr()];
    }


    public void redrawAncestralStates() {
        jointLogLikelihood = 0;
        TreeInterface tree = treeInput.get();
        traverseSample(tree, tree.getRoot(), null);

        areStatesRedrawn = true;
    }

//    private boolean checkConditioning = true;


    @Override
    public double calculateLogP() {
        areStatesRedrawn = false;
        double marginalLogLikelihood = super.calculateLogP();
        likelihoodKnown = true;

        if (returnMarginalLogLikelihood) {
            return marginalLogLikelihood;
        }
        // redraw states and return joint density of drawn states
        redrawAncestralStates();
        logP = jointLogLikelihood;
        return logP;
    }

    protected Helper treeTraits = new Helper();

    public TreeTrait[] getTreeTraits() {
        return treeTraits.getTreeTraits();
    }

    public TreeTrait getTreeTrait(String key) {
        return treeTraits.getTreeTrait(key);
    }


    private static String formattedState(int[] state, DataType dataType) {
        StringBuilder sb = new StringBuilder();
        sb.append("\"");
        if (dataType instanceof UserDataType) {
            boolean first = true;
            for (int i : state) {
                if (!first) {
                    sb.append(" ");
                } else {
                    first = false;
                }

                sb.append(dataType.getCode(i));
            }

        } else {
            for (int i : state) {
                sb.append(dataType.getChar(i));
            }
        }
        sb.append("\"");
        return sb.toString();
    }

    private int drawChoice(double[] measure) {
        if (useMAP) {
            double max = measure[0];
            int choice = 0;
            for (int i = 1; i < measure.length; i++) {
                if (measure[i] > max) {
                    max = measure[i];
                    choice = i;
                }
            }
            return choice;
        } else {
            return Randomizer.randomChoicePDF(measure);
        }
    }

    public void getStates(int tipNum, int[] states)  {
        // Saved locally to reduce BEAGLE library access
        System.arraycopy(tipStates[tipNum], 0, states, 0, states.length);
    }


    /**
     * Traverse (pre-order) the tree sampling the internal node states.
     *
     * @param tree        - TreeModel on which to perform sampling
     * @param node        - current node
     * @param parentState - character state of the parent node to 'node'
     */
    public void traverseSample(TreeInterface tree, Node node, int[] parentState) {

        int nodeNum = node.getNr();

        Node parent = node.getParent();

        // This function assumes that all partial likelihoods have already been calculated
        // If the node is internal, then sample its state given the state of its parent (pre-order traversal).

        double[] conditionalProbabilities = new double[stateCount];
        int[] state = new int[patternCount];

        if (!node.isLeaf()) {

            if (parent == null) {

                double[] rootPartials = m_fRootPartials;

                double[] rootFrequencies = substitutionModel.getFrequencies();
                if (rootFrequenciesInput.get() != null) {
                    rootFrequencies = rootFrequenciesInput.get().getFreqs();
                }

                // This is the root node
                for (int j = 0; j < patternCount; j++) {
                	if (beagle != null) {
                        System.err.println("beagle support not implemented yet");
                        System.exit(1);
                	} else {
                		System.arraycopy(rootPartials, j * stateCount, conditionalProbabilities, 0, stateCount);
                	}

                    for (int i = 0; i < stateCount; i++) {
                        conditionalProbabilities[i] *= rootFrequencies[i];
                    }
                    try {
                        state[j] = drawChoice(conditionalProbabilities);
                    } catch (Error e) {
                        System.err.println(e);
                        System.err.println("Please report error to Marc");
                        state[j] = 0;
                    }
                    reconstructedStates[nodeNum][j] = state[j];
                    jointLogLikelihood += Math.log(rootFrequencies[state[j]]);
                }

            } else {

                // This is an internal node, but not the root
                double[] partialLikelihood = new double[stateCount * patternCount];

            	if (beagle != null) {
                    System.err.println("beagle support not implemented yet");
                    System.exit(1);
            	} else {
                    likelihoodCore.getNodePartials(node.getNr(), partialLikelihood);
                    likelihoodCore.getNodeMatrix(nodeNum, 0, probabilities);
            	}


                for (int j = 0; j < patternCount; j++) {

                    int parentIndex = parentState[j] * stateCount;
                    int childIndex = j * stateCount;

                    for (int i = 0; i < stateCount; i++) {
                        conditionalProbabilities[i] = partialLikelihood[childIndex + i] * probabilities[parentIndex + i];
                    }

                    state[j] = drawChoice(conditionalProbabilities);
                    reconstructedStates[nodeNum][j] = state[j];
                    double contrib = probabilities[parentIndex + state[j]];
                    jointLogLikelihood += Math.log(contrib);
                }
            }

            // Traverse down the two child nodes
            Node child1 = node.getChild(0);
            traverseSample(tree, child1, state);

            Node child2 = node.getChild(1);
            traverseSample(tree, child2, state);
        } else {

            // This is an external leaf
        	getStates(nodeNum, reconstructedStates[nodeNum]);

        	if (sampleTipsInput.get()) {
	            // Check for ambiguity codes and sample them
	            for (int j = 0; j < patternCount; j++) {

	                final int thisState = reconstructedStates[nodeNum][j];
	                final int parentIndex = parentState[j] * stateCount;
	            	if (beagle != null) {
                        System.err.println("beagle support not implemented yet");
                        System.exit(1);
	            	} else {
	                likelihoodCore.getNodeMatrix(nodeNum, 0, probabilities);
	            	}
	                if (dataType.isAmbiguousCode(thisState)) {

	                    boolean [] stateSet = dataType.getStateSet(thisState);
	                    for (int i = 0; i < stateCount; i++) {
	                        conditionalProbabilities[i] =  stateSet[i] ? probabilities[parentIndex + i] : 0;
	                    }

	                    reconstructedStates[nodeNum][j] = drawChoice(conditionalProbabilities);
	                }

	                double contrib = probabilities[parentIndex + reconstructedStates[nodeNum][j]];
	                jointLogLikelihood += Math.log(contrib);
	            }
        	}

        }
    }


    @Override
    public void log(final long sample, final PrintStream out) {
    	// useful when logging on a fixed tree in an AncestralTreeLikelihood that is logged, but not part of the posterior
    	hasDirt = Tree.IS_FILTHY;
    	calculateLogP();
        out.print(getCurrentLogP() + "\t");
    }

    @Override
    public void finalize() throws Throwable {
        getLikelihoodCore().finalize();
        reconstructedStates = null;
        storedReconstructedStates = null;
        probabilities = null;
        tag = null;
        jointLogLikelihood = 0;
        storedJointLogLikelihood = 0;
        patternCount=0;
        stateCount = 0;
        m_fRootPartials = null;
        substitutionModel = null;
        tipStates = null;
        logP = 0;
        traitDimension = 0;
        leafNr = null;
        storedTipStates = null;
        parameters = null;
        treeTraits = null;
        treeTraits = new Helper();
    }
}