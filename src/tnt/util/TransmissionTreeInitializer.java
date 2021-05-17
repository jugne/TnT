package tnt.util;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import beast.core.BEASTInterface;
import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.alignment.distance.Distance;
import beast.evolution.alignment.distance.JukesCantorDistance;
import beast.evolution.tree.Node;
import beast.evolution.tree.RandomTree;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.math.distributions.MRCAPrior;
import beast.util.ClusterTree;
import beast.util.Randomizer;
import beast.util.TreeParser;
import starbeast2.PopulationModel;
import starbeast2.SpeciesTreeInterface;
import tnt.transmissionTree.TransmissionTree;

/**
 * Adapted by Ugne Stolz, to keep tree direction correct when parsing from
 * pre-defined newick string and allow for same patient sampling over time.
 * 
 * @author Joseph Heled
 * @author Huw Ogilvie
 */

@Description("Set a starting point for a *BEAST analysis from gene alignment data.")
public class TransmissionTreeInitializer extends Tree implements StateNodeInitialiser {

    static enum Method {
        POINT("point-estimate"),
        ALL_RANDOM("random");

        Method(final String name) {
            this.ename = name;
        }

        @Override
		public String toString() {
            return ename;
        }

        private final String ename;
    }
    final public Input<Method> initMethod = new Input<>("method", "Initialise either with a totally random " +
            "state or a point estimate based on alignments data (default point-estimate)",
            Method.POINT, Method.values());

	final public Input<TransmissionTree> transmissionTreeInput = new Input<>("transmissionTree",
			"The species tree to initialize.",
			Input.Validate.REQUIRED);
    final public Input<String> newickInput = new Input<>("newick", "Newick string for a custom initial species tree.");

	public Input<List<TaxonSet>> taxonsetsInput = new Input<>("patientTaxonSets",
			"a separate list of taxa for samples collected from the same patient", new ArrayList<>());

//    final public Input<List<Tree>> genes = new Input<>("geneTree", "Gene trees to initialize", new ArrayList<>());

	final public Input<List<GeneTreeInitializer>> genesTreeInitializerInput = new Input<>("geneTreeInitializerList",
			"New gene tree initializer for each tree", new ArrayList<>());

    final public Input<RealParameter> birthRate = new Input<>("birthRate",
            "Tree prior birth rate to initialize");

    final public Input<Function> muInput = new Input<>("baseRate",
            "Main clock rate used to scale trees (default 1).");
    
//	public Input<TraitSet> sampleCountsInput = new Input<>(
//			"sampleCounts",
//			"TraitSet defining number of  samples per node in species tree.");

	public Input<RealParameter> popSizesInput = new Input<RealParameter>("populationSizes",
			"Constant per-branch effective population sizes.");

	public Input<RealParameter> originInput = new Input<RealParameter>("origin", "origin");
	public Input<Function> deathRateInput = new Input<Function>("deathRate", "Death rate");
	public Input<RealParameter> samplingRateInput = new Input<RealParameter>("samplingRate",
			"Sampling rate per individual");

	public Input<RealParameter> bottleneckStrengthInput = new Input<RealParameter>("bottleneckStrength",
			"Strength of the bottleneck in scaled time");


	final public Input<PopulationModel> populationFunctionInput = new Input<>("populationModel",
			"The species tree population model.");

	private Map<String, String> tipSpeciesMap;
	private Map<String, Node> speciesNodeMap;
	private Set<String> allSpeciesNames;
	private Set<String> allTipNames;


    @Override
    public void initStateNodes() {
		Log.info.println(
				"TnT transmissionTree initialiser credits StarBEAST2 speciesTree initializer for most of its functionality!");

		final TransmissionTree transmissionTree = transmissionTreeInput.get();
		final TaxonSet taxonSuperSet = transmissionTree.getTaxonset();
		final Set<BEASTInterface> treeOutputs = transmissionTreeInput.get().getOutputs();
        final Method method = initMethod.get();
        final String newick = newickInput.get();

		tipSpeciesMap = new HashMap<>();
		speciesNodeMap = new HashMap<>();
		allTipNames = new HashSet<>();
		allSpeciesNames = new HashSet<>();

		for (Taxon species : taxonSuperSet.taxonsetInput.get()) {
			final String speciesName = species.getID();
			final TaxonSet speciesTaxonSet = (TaxonSet) species;

			allSpeciesNames.add(speciesName);

			for (Taxon tip : speciesTaxonSet.taxonsetInput.get()) {
				final String tipName = tip.getID();
				tipSpeciesMap.put(tipName, speciesName);
				allTipNames.add(tipName);
			}
		}

		for (Node node : transmissionTree.getExternalNodes())
			speciesNodeMap.put(node.getID(), node);

        final List<MRCAPrior> calibrations = new ArrayList<>();
        for (final Object plugin : treeOutputs ) {
            if (plugin instanceof MRCAPrior) {
                calibrations.add((MRCAPrior) plugin);
            }
        }

		final List<GeneTreeInitializer> geneTrees = genesTreeInitializerInput.get();
        boolean userSpecifiedGeneTrees = false;
		for (final GeneTreeInitializer gtree : geneTrees) {
			if (gtree.getGeneTree() instanceof beast.util.TreeParser) {
                userSpecifiedGeneTrees = true;
                break;
            }
        }



        boolean geneTreesNeedInit = true;
        if (newick != null) {
			Log.info.println("TnT: using newick string to initialize transmission tree.");
			final TreeParser parser = new TreeParser(newick);
			copyTreeStructure(parser, transmissionTree);
		} else if (!taxonsetsInput.get().isEmpty()) {
			Log.info.println("TnT: using same patient sampling to initialize transmission tree.");
			samePatientSamplingInit(transmissionTree);
        } else if (method == Method.ALL_RANDOM) {
			Log.info.println("TnT: using randomInit to initialize transmission tree.");
			randomInit(transmissionTree, calibrations);
        } else if (!checkSpeciesAlwaysRepresented())  {
			Log.info.println("TnT: using randomInit to initialize transmission tree (required by missing data)).");
			randomInit(transmissionTree, calibrations);
        } else if (calibrations.size() > 0)  {
			Log.info.println("TnT: using randomInit to initialize transmission tree (required by calibrations)).");
			randomInit(transmissionTree, calibrations);
		} else if (transmissionTree.hasDateTrait()) {
			Log.info.println("TnT: using randomInit to initialize transmission tree (required by tip dates).");
			randomInit(transmissionTree, calibrations);
        } else if (userSpecifiedGeneTrees) {
			Log.info.println(
					"TnT: using randomInit to initialize transmission tree (required by user-specified gene trees).");
			randomInit(transmissionTree, calibrations);
        } else if (method == Method.POINT) {
			Log.info.println("TnT: using fullInit to initialize all trees.");
			fullInit(transmissionTree);
            geneTreesNeedInit = false;
        }

        if (geneTreesNeedInit) {
			final double rootHeight = transmissionTree.getRoot().getHeight();
			Log.info.println(String
					.format("TnT: initializing gene trees with random or user-specified topologies (%f).", rootHeight));

			for (final GeneTreeInitializer gtree : geneTrees) {
				if (gtree.getGeneTree() instanceof beast.util.TreeParser) {
                    /* add the height of the species root node to all gene tree node height, to
                    ensure compatibility of the trees while preserving user-specified topologies */
					boostGeneTreeInternalNodeHeights(gtree.getGeneTree(), rootHeight);
				} else if (bottleneckStrengthInput.get() != null && gtree.getSampleCounts() != null) {
					final tnt.simulator.SimulatedGeneTreeInit geneTree = new tnt.simulator.SimulatedGeneTreeInit();
					gtree.getGeneTree().setID("gene_tree_truth");
					gtree.getGeneTree().initByName("transmissionTreeInput", transmissionTree,
							"sampleCounts", gtree.getSampleCounts(),
							"popSizeAboveOrigin",
							new RealParameter(Double.toString(popSizesInput.get().getArrayValue(1))),
							"populationSizes", new RealParameter(Double.toString(popSizesInput.get().getArrayValue(0))),
							"bottleneckStrength", bottleneckStrengthInput.get(),
							"birthRate", birthRate.get(),
							"deathRate", deathRateInput.get(),
							"samplingRate", samplingRateInput.get(),
							"origin", originInput.get());

//					final TreeParser parser = new TreeParser(geneTree.getSimulatedGeneTree().getRoot().toNewick());
//					gtree.assignFromWithoutID(parser);
					gtree.getGeneTree().assignFromWithoutID(geneTree.getSimulatedGeneTree());
				} else {
					gtree.getGeneTree().makeCaterpillar(rootHeight,
							rootHeight / gtree.getGeneTree().getInternalNodeCount(), true);
                }
				// make sure the heights of all gene tree tips is equal to the height of
				// corresponding species tree tips
				resetGeneTreeTipHeights(transmissionTree, gtree.getGeneTree());
            }
        }

        // initialize population sizes to equal average branch length
        // this is equivalent to 2Ne = E[1/lambda]
		final double speciesTreeLength = TreeStats.getLength(transmissionTree);
		final int nBranches = transmissionTree.getNodeCount();
        final double averageBranchLength = speciesTreeLength / (nBranches - 1);

        final PopulationModel populationModel = populationFunctionInput.get();
        if (populationModel != null) populationModel.initPopSizes(averageBranchLength);
        
		for (GeneTreeInitializer geneTree : geneTrees) {
			for (Node geneNode : geneTree.getGeneTree().getNodesAsArray()) {
        		if (geneNode.isLeaf()) {
        			final String tipName = geneNode.getID();
					if (!allTipNames.contains(tipName)) {
        	            throw new RuntimeException(String.format("ERROR: Gene tree tip name '%s' is missing from taxon map. "
        	            		+ "This typically occurs when a sequence or sample name is identical to a species name. "
        	            		+ "Make sure all species names are distinct from sequence or sample names.", tipName));
        			}
        		}
        	}
        }
    }

    private void boostGeneTreeInternalNodeHeights(Tree gtree, double boost) {
        for (Node node: gtree.getNodesAsArray()) {
            if (!node.isLeaf()) {
                final double newHeight = node.getHeight() + boost;
                node.setHeight(newHeight);
            }
        }
    }

    private boolean checkSpeciesAlwaysRepresented() {
		for (GeneTreeInitializer geneTree : genesTreeInitializerInput.get()) {
			final String[] allTipNames = geneTree.getGeneTree().getTaxaNames();

			for (String speciesName : allSpeciesNames) {
				boolean speciesRepresented = false;
				for (String tipName : allTipNames)
					speciesRepresented |= tipSpeciesMap.get(tipName).equals(speciesName);
				if (!speciesRepresented)
					return false;
			}
        }

        return true;
    }

	private void resetGeneTreeTipHeights(SpeciesTreeInterface speciesTree, Tree gtree) {
		for (Node geneLeaf : gtree.getExternalNodes()) {
			if (allTipNames.contains(geneLeaf.getID())) {
				final String speciesName = tipSpeciesMap.get(geneLeaf.getID());
				final Node speciesLeaf = speciesNodeMap.get(speciesName);
				geneLeaf.setHeight(speciesLeaf.getHeight());
			} else {
				throw new RuntimeException("The taxon " + geneLeaf.getID() + " is missing from the taxonsuperset!");
			}
        }
	}

	private double[] firstMeetings(final Tree gtree, final Map<String, Integer> tipName2Species, final int speciesCount) {
        final Node[] nodes = gtree.listNodesPostOrder(null, null);
        @SuppressWarnings("unchecked")
		final Set<Integer>[] tipsSpecies = new Set[nodes.length];
        for(int k = 0; k < tipsSpecies.length; ++k) {
            tipsSpecies[k] = new LinkedHashSet<>();
        }
        // d[i,j] = minimum height of node which has tips belonging to species i and j
        // d is is upper triangular
        final double[] dmin = new double[(speciesCount*(speciesCount-1))/2];
        Arrays.fill(dmin, Double.MAX_VALUE);

        for (final Node n : nodes) {
            if (n.isLeaf()) {
                tipsSpecies[n.getNr()].add(tipName2Species.get(n.getID()));
            } else {
                assert n.getChildCount() == 2;
                @SuppressWarnings("unchecked")
				final Set<Integer>[] sps = new Set[2];
                sps[0] = tipsSpecies[n.getChild(0).getNr()];
                sps[1] = tipsSpecies[n.getChild(1).getNr()];
                final Set<Integer> u = new LinkedHashSet<>(sps[0]);
                u.retainAll(sps[1]);
                sps[0].removeAll(u);
                sps[1].removeAll(u);

                for (final Integer s1 : sps[0]) {
                    for (final Integer s2 : sps[1]) {
                        final int i = getDMindex(speciesCount, s1, s2);
                        dmin[i] = min(dmin[i], n.getHeight());
                    }
                }
                u.addAll(sps[0]);
                u.addAll(sps[1]);
                tipsSpecies[n.getNr()] = u;
            }
        }
        return dmin;
    }

    private int getDMindex(final int speciesCount, final int s1, final int s2) {
        final int mij = min(s1,s2);
        return (mij*(2*speciesCount-1 - mij))/2 + (abs(s1-s2)-1);
    }


	private void fullInit(final TransmissionTree transmissionTree) {
        // Build gene trees from  alignments
    	

        final Function muInput = this.muInput.get();
        final double mu =  (muInput != null )  ? muInput.getArrayValue() : 1;

        final TaxonSet species = transmissionTree.m_taxonset.get();
        final List<String> speciesNames = species.asStringList();
        final int speciesCount = speciesNames.size();

		final List<GeneTreeInitializer> geneTrees = genesTreeInitializerInput.get();

        //final List<Alignment> alignments = genes.get();
        //final List<Tree> geneTrees = new ArrayList<>(alignments.size());
        double maxNsites = 0;
        //for( final Alignment alignment : alignments)  {
		for (final GeneTreeInitializer gtree : geneTrees) {
            //final Tree gtree = new Tree();
			final Alignment alignment = gtree.getGeneTree().m_taxonset.get().alignmentInput.get();

            final ClusterTree ctree = new ClusterTree();
			ctree.initByName("initial", gtree.getGeneTree(), "clusterType", "upgma", "taxa", alignment);
			gtree.getGeneTree().scale(1 / mu);

            maxNsites = max(maxNsites, alignment.getSiteCount());
        }
        final Map<String, Integer> geneTips2Species = new LinkedHashMap<>();
        final List<Taxon> taxonSets = species.taxonsetInput.get();

        for(int k = 0; k < speciesNames.size(); ++k) {
            final Taxon nx = taxonSets.get(k);
            final List<Taxon> taxa = ((TaxonSet) nx).taxonsetInput.get();
            for( final Taxon n : taxa ) {
              geneTips2Species.put(n.getID(), k);
            }
        }
        final double[] dg = new double[(speciesCount*(speciesCount-1))/2];

        final double[][] genesDmins = new double[geneTrees.size()][];

        for( int ng = 0; ng < geneTrees.size(); ++ng ) {
			final Tree g = geneTrees.get(ng).getGeneTree();
            final double[] dmin = firstMeetings(g, geneTips2Species, speciesCount);
            genesDmins[ng] = dmin;

            for(int i = 0; i < dmin.length; ++i) {
                dg[i] += dmin[i];
                if (dmin[i] == Double.MAX_VALUE) {
                	// this happens when a gene tree has no taxa for some species-tree taxon.
                	// TODO: ensure that if this happens, there will always be an "infinite"
                	// distance between species-taxon 0 and the species-taxon with missing lineages,
                	// so i < speciesCount - 1.
                	// What if lineages for species-taxon 0 are missing? Then all entries will be 'infinite'.
                	String id = (i < speciesCount - 1? transmissionTree.getExternalNodes().get(i+1).getID() : "unknown taxon");
                	if (i == 0) {
                		// test that all entries are 'infinite', which implies taxon 0 has lineages missing 
                		boolean b = true;
                		for (int k = 1; b && k < speciesCount - 1; k++) {
                			b = (dmin[k] == Double.MAX_VALUE);
                		}
                		if (b) {
                			// if all entries have 'infinite' distances, it is probably the first taxon that is at fault
                			id = transmissionTree.getExternalNodes().get(0).getID();
                		}
                	}
                	throw new RuntimeException("Gene tree " + g.getID() + " has no lineages for species taxon " + id + " ");
                }
            }
        }

        for(int i = 0; i < dg.length; ++i) {
            double d = dg[i] / geneTrees.size();
            if( d == 0 ) {
               d = (0.5/maxNsites) * (1/mu);
            } else {
                // heights to distances
                d *= 2;
            }
            dg[i] = d;
        }

        final ClusterTree ctree = new ClusterTree();
        final Distance distance = new Distance() {
            @Override
            public double pairwiseDistance(final int s1, final int s2) {
                final int i = getDMindex(speciesCount, s1,s2);
                return dg[i];
            }
        };
        ctree.initByName("initial", transmissionTree, "taxonset", species,"clusterType", "upgma", "distance", distance);

        final Map<String, Integer> sptips2SpeciesIndex = new LinkedHashMap<>();
        for(int i = 0; i < speciesNames.size(); ++i) {
            sptips2SpeciesIndex.put(speciesNames.get(i), i);
        }
        final double[] spmin = firstMeetings(transmissionTree, sptips2SpeciesIndex, speciesCount);

        for( int ng = 0; ng < geneTrees.size(); ++ng ) {
            final double[] dmin = genesDmins[ng];
            boolean compatible = true;
            for(int i = 0; i < spmin.length; ++i) {
                if( dmin[i] <= spmin[i] ) {
                    compatible = false;
                    break;
                }
            }
            if( ! compatible ) {
				final Tree gtree = geneTrees.get(ng).getGeneTree();
                final TaxonSet gtreeTaxa = gtree.m_taxonset.get();
                final Alignment alignment = gtreeTaxa.alignmentInput.get();
                final List<String> taxaNames = alignment.getTaxaNames();
                final int taxonCount =  taxaNames.size();
                // speedup
                final Map<Integer,Integer> g2s = new LinkedHashMap<>();
                for(int i = 0; i < taxonCount; ++i) {
                    g2s.put(i, geneTips2Species.get(taxaNames.get(i)));
                }

                final JukesCantorDistance jc = new JukesCantorDistance();
                jc.setPatterns(alignment);
                final Distance gdistance = new Distance() {
                    @Override
                    public double pairwiseDistance(final int t1, final int t2) {
                        final int s1 = g2s.get(t1);
                        final int s2 = g2s.get(t2);
                        double d = jc.pairwiseDistance(t1,t2)/mu;
                        if( s1 != s2 ) {
                            final int i = getDMindex(speciesCount, s1,s2);
                            final double minDist = 2 * spmin[i];
                            if( d <= minDist ) {
                                d = minDist * 1.001;
                            }
                        }
                        return d;
                    }
                };
                final ClusterTree gtreec = new ClusterTree();
                gtreec.initByName("initial", gtree, "taxonset", gtreeTaxa,
                        "clusterType", "upgma", "distance", gdistance);
            }
        }

        final RealParameter lambda = birthRate.get();
        if (lambda != null && lambda instanceof StateNode) {
            final StateNode lambdaStateNode = lambda;

            // only change lambda if it is to be estimated
            if (lambdaStateNode.isEstimatedInput.get()) {
                final double rh = transmissionTree.getRoot().getHeight();
                double l = 0;
                for(int i = 2; i < speciesCount+1; ++i) l += 1./i;
                lambda.setValue((1 / rh) * l);
            }
        }
    }

	private void randomInit(final TransmissionTree transmissionTree, List<MRCAPrior> calibrations) {
    	final RealParameter birthRateParameter = birthRate.get();
    	final Double lambda = (birthRateParameter == null) ? 1.0 : birthRateParameter.getValue();
    	final Double initialPopSize = 1.0 / lambda; // scales coalescent tree height inverse to birth rate
    	final RealParameter popSize = new RealParameter(initialPopSize.toString());
        final ConstantPopulation pf = new ConstantPopulation();
        pf.setInputValue("popSize", popSize);

        final RandomTree rnd = new RandomTree();
        rnd.setInputValue("taxonset", transmissionTree.getTaxonset());
        if (transmissionTree.hasDateTrait()) rnd.setInputValue("trait", transmissionTree.getDateTrait());

        for (final MRCAPrior cal: calibrations) rnd.setInputValue("constraint", cal);

        rnd.setInputValue("populationModel", pf);
        rnd.setInputValue("populationModel", pf);
        rnd.initAndValidate();

		copyTreeStructure(rnd, transmissionTree);
	}

	private void samePatientSamplingInit(TransmissionTree transmissionTree) {
		final RealParameter birthRateParameter = birthRate.get();
		final Double lambda = (birthRateParameter == null) ? 1.0 : birthRateParameter.getValue();
		final Double initialPopSize = 1.0 / lambda; // scales coalescent tree height inverse to birth rate
		final RealParameter popSize = new RealParameter(initialPopSize.toString());
		final ConstantPopulation pf = new ConstantPopulation();
		pf.setInputValue("popSize", popSize);

		final RandomTree rnd = new RandomTree();
		rnd.setInputValue("taxonset", transmissionTree.getTaxonset());
		if (transmissionTree.hasDateTrait())
			rnd.setInputValue("trait", transmissionTree.getDateTrait());

		rnd.setInputValue("populationModel", pf);
		rnd.setInputValue("populationModel", pf);
		rnd.initAndValidate();

		System.out.println(rnd.getRoot().toNewick());

		/////////////

		int nSets = taxonsetsInput.get().size();
		int nSamples = transmissionTree.getLeafNodeCount();
		List<List<Node>> taxonsets = new ArrayList<List<Node>>(nSets);
		List<Node> singlePatientSamples = new ArrayList<Node>();
		double maxHeight = Double.NEGATIVE_INFINITY;

		// If there are same patient sampling over time supplied
		if (nSets != 0) {
			for (int i = 0; i < nSets; i++)
				taxonsets.add(new ArrayList<Node>());

			for (int idx = 0; idx < nSamples; idx++) {
				Node leaf = rnd.getNode(idx);
				if (leaf.getParent() != null)
					leaf.getParent().removeAllChildren(false);
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
				singlePatientSamples.add(leaf);
				if (leaf.getHeight() > maxHeight)
					maxHeight = leaf.getHeight();
			}

			for (int ii = 0; ii < nSets; ii++) {
				taxonsets.get(ii).sort(Comparator.comparing(Node::getHeight));
			}
		}

		int id = nSamples;
		for (int s = 0; s < nSets; s++) {
			Node[] tmp = new Node[taxonsets.get(s).size()];
			tmp = taxonsets.get(s).toArray(tmp);
			int nNodes = tmp.length;
			for (int n = 0; n < nNodes - 1; n++) {
				Node parent = new Node();
				parent.setID(Integer.toString(id));
				parent.setNr(id);
				id += 1;

				parent.addChild(tmp[n]);
				parent.addChild(tmp[n + 1]);
				singlePatientSamples.remove(tmp[n]);
				singlePatientSamples.remove(tmp[n + 1]);
				parent.setHeight(tmp[n + 1].getHeight());
				tmp[n] = null;
				tmp[n + 1] = parent;
			}
			singlePatientSamples.add(tmp[nNodes - 1]);
			if (tmp[nNodes - 1].getHeight() > maxHeight)
				maxHeight = tmp[nNodes - 1].getHeight();
		}

		while (singlePatientSamples.size() > 1) {
			Node node1 = singlePatientSamples.get(Randomizer.nextInt(singlePatientSamples.size()));
			singlePatientSamples.remove(node1);
			Node node2 = singlePatientSamples.get(Randomizer.nextInt(singlePatientSamples.size()));
			singlePatientSamples.remove(node2);

			double deltaT = Randomizer.nextExponential(
					(singlePatientSamples.size() + 2) * (singlePatientSamples.size() + 1) * 0.5 * 1.0 / initialPopSize);
			double coalTime = maxHeight + deltaT;

			Node parent = new Node();
			parent.setID(Integer.toString(id));
			parent.setNr(id);
			id += 1;
			parent.setHeight(coalTime);
			parent.addChild(node1);
			parent.addChild(node2);

			singlePatientSamples.add(parent);
			maxHeight = coalTime;
		}

		Node root = singlePatientSamples.get(0);

		copyTreeStructure(new Tree(root), transmissionTree);
	}

	// copy the structure of the source tree to the destination tree
	// preserving the leaf node names and numbers
	private void copyTreeStructure(final Tree src, final Tree dst) {
		final Node[] dstNodes = dst.getNodesAsArray();
		final Node[] srcNodes = src.getNodesAsArray();

		final Map<String, Integer> srcTipNumbers = new HashMap<>();

		final int nodeCount = src.getNodeCount();
		final int leafNodeCount = src.getLeafNodeCount();

		// Clear the children of all internal nodes in the destination tree
		for (int nodeNumber = leafNodeCount; nodeNumber < nodeCount; nodeNumber++)
			dstNodes[nodeNumber].removeAllChildren(false);

		// Record the node number of all leaves in the source tree
		for (int nodeNumber = 0; nodeNumber < leafNodeCount; nodeNumber++) {
			final String srcName = srcNodes[nodeNumber].getID();
			srcTipNumbers.put(srcName, nodeNumber);
		}

		// Set the heights of all nodes to match the source height
		for (int nodeNumber = 0; nodeNumber < nodeCount; nodeNumber++) {
			final Node dstNode = dstNodes[nodeNumber];
			Node srcNode;

			// find the corresponding node from the source tree
			if (nodeNumber < leafNodeCount) { // if this is a leaf node
				final String speciesName = dstNode.getID();
				System.out.println(speciesName);
				final int srcTipNumber = srcTipNumbers.get(speciesName);

				srcNode = srcNodes[srcTipNumber];
			} else { // if this is an internal node
				srcNode = srcNodes[nodeNumber];
			}

			// Copy height
			dstNode.setHeight(srcNode.getHeight());

			// Clear and copy metadata
			final Set<String> dstMetaDataNames = dstNode.getMetaDataNames();
			final Set<String> srcMetaDataNames = srcNode.getMetaDataNames();
			final Set<String> srcLengthMetaDataNames = srcNode.getLengthMetaDataNames();

			for (String metaDataName : dstMetaDataNames)
				dstNode.removeMetaData(metaDataName);

			for (String metaDataName : srcMetaDataNames)
				dstNode.setMetaData(metaDataName, srcNode.getMetaData(metaDataName));

			for (String lengthMetaDataName : srcLengthMetaDataNames)
				dstNode.setMetaData(lengthMetaDataName, srcNode.getLengthMetaData(lengthMetaDataName));

			// if this is not the root node, also set the parent and child
			// connections to match the source
			if (nodeNumber != nodeCount - 1) {
				boolean rightChild = false;
				Node siblingNode = null;
				final int parentNumber = srcNode.getParent().getNr();
				final Node parentNode = dstNodes[parentNumber];
				if (srcNode.getParent().getChild(0) == srcNode)
					rightChild = true;
				if (rightChild && parentNode.getChildCount() > 0) {
					siblingNode = parentNode.getChild(0);
					parentNode.removeChild(siblingNode);
				}
				dstNode.setParent(parentNode);
				parentNode.addChild(dstNode);

				if (siblingNode != null) {
					siblingNode.setParent(parentNode);
					parentNode.addChild(siblingNode);
				}

			}
		}
    }

    @Override
    public void getInitialisedStateNodes(final List<StateNode> stateNodes) {
		stateNodes.add(transmissionTreeInput.get());

		for (final GeneTreeInitializer g : genesTreeInitializerInput.get()) {
			stateNodes.add(g.getGeneTree());
        }

        final RealParameter brate = birthRate.get();
        if (brate != null) {
            stateNodes.add(brate) ;
        }
    }
}
