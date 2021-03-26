package tnt.orientation;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.junit.Test;

import beast.core.BEASTObject;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Logger;
import beast.core.MCMC;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.core.util.CompoundDistribution;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.speciation.SABirthDeathModel;
import beast.evolution.substitutionmodel.JukesCantor;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.math.distributions.Prior;
import beast.math.distributions.Uniform;
import beast.util.Randomizer;
import junit.framework.Assert;
//import starbeast2.NodeReheight2;
import starbeast2.SpeciesTreeInterface;
import starbeast2.StarBeastTaxonSet;
import tnt.distribution.GeneTreeDistribution;
import tnt.operators.CoordinatedExponential;
import tnt.operators.CoordinatedUniform;
import tnt.operators.CreateMergersOrReheight;
import tnt.operators.ExchangeOperator;
import tnt.operators.ExpandCollapseOperator;
//import beast.evolution.operators.Exchange;
import tnt.operators.LeafToSampledAncestorJump;
import tnt.operators.SAUniform;
import tnt.operators.SAWilsonBalding;
//import beast.evolution.operators.SAExchange;
//import beast.evolution.operators.SAScaleOperator;
import tnt.operators.SPROperatorChanged;
import tnt.operators.TransmissionAttachChanged;

public class OrientationTest{
	
	String initNewick = "(1:0.05244677605979664,(2:0.21990168652563447,3:1.2199016865256345):0.8325450895341622):0.0;";

	@Test
	public void topologyDistribution() throws Exception {
		
		Randomizer.setSeed(127);
		
		// data
		List<Object> alignmentInitArgs = new ArrayList<Object>();
		for (int i = 1; i <= 3; i++) {
			Sequence thisSeq = new Sequence();
			thisSeq.initByName("taxon", String.valueOf(i) + "_x", "totalcount",
					"4", "value", "-");
			thisSeq.setID("seq_" + String.valueOf(i) + "_x");
			alignmentInitArgs.add("sequence");
			alignmentInitArgs.add(thisSeq);
		}
		Alignment alignment = new Alignment();
		alignment.initByName(alignmentInitArgs.toArray());
		alignment.setID("Gene");

		// state

		TaxonSet taxonSet1 = new TaxonSet();
		Taxon taxonSet11 = new Taxon();
		taxonSet1.setID("1");
		taxonSet11.setID("1_x");
		taxonSet1.initByName("taxon", taxonSet11);

		TaxonSet taxonSet2 = new TaxonSet();
		Taxon taxonSet21 = new Taxon();
		taxonSet2.setID("2");
		taxonSet21.setID("2_x");
		taxonSet2.initByName("taxon", taxonSet21);

		TaxonSet taxonSet3 = new TaxonSet();
		Taxon taxonSet31 = new Taxon();
		taxonSet3.setID("3");
		taxonSet31.setID("3_x");
		taxonSet3.initByName("taxon", taxonSet31);

		StarBeastTaxonSet starTaxonSet = new StarBeastTaxonSet();
		tnt.transmissionTree.TransmissionTree transmissionTree = new tnt.transmissionTree.TransmissionTree();

		starTaxonSet.initByName("taxon", taxonSet1, "taxon", taxonSet2, "taxon", taxonSet3);
		starTaxonSet.setID("taxonsuperset");

		TraitSet trait = new TraitSet();
		trait.initByName("traitname", "date-backward", "value", "1=2,2=1,3=0", "taxa", starTaxonSet);
		trait.setID("dateTrait.t:Species");

		transmissionTree.initByName("trait", trait, "taxonset", starTaxonSet);
		transmissionTree.setID("Tree.t:Species");

		TaxonSet taxonSet = new TaxonSet();
		taxonSet.initByName("alignment", alignment);

		Tree geneTree = new Tree();
		geneTree.setID("Tree.t:Gene");
		geneTree.initByName("taxonset", taxonSet);

		RealParameter origin = new RealParameter("3.0");
		origin.initByName("lower", "0.0", "upper", "Infinity");
		origin.setID("originFBD.t:Species");

		// Set up state:
		State state = new State();
		state.initByName("stateNode", transmissionTree, "stateNode", geneTree, "stateNode", origin);
		state.setID("state");

		// Set up population size and bottleneck:
		RealParameter popSize = new RealParameter("0.5");
		popSize.initByName("lower", 0.0, "upper", 2.0, "estimate", false);
		popSize.setID("constPopSizes.Species");

		RealParameter tau = new RealParameter("0.5");
		popSize.initByName("lower", 0.1, "upper", 1.0, "estimate", false);
		popSize.setID("constPopSizes.Species");

		
		tnt.distribution.GeneTreeIntervals intervals = new tnt.distribution.GeneTreeIntervals();
		intervals.initByName("geneTree", geneTree, "transmissionTree", transmissionTree);

		// Set up distributions:

		// dist 1
		CompoundDistribution prior = new CompoundDistribution();
		SABirthDeathModel FBDModel = new SABirthDeathModel();
		Prior priorDist = new Prior();
		Uniform uniform3 = new Uniform();


		RealParameter birthRate = new RealParameter("2.0");
		RealParameter deathRate = new RealParameter("1.0");
		RealParameter samplingRate = new RealParameter("0.5");
		RealParameter removalProbability = new RealParameter("0.9");

		birthRate.initByName("estimate", false, "lower", "0.0");
		deathRate.initByName("estimate", false, "lower", "0.0", "upper", "1.0");
		samplingRate.initByName("estimate", false, "lower", "0.0", "upper", "1.0");
		removalProbability.initByName("lower", "0.0", "upper", "1.0");

		starbeast2.StarBeastInitializer initializer = new starbeast2.StarBeastInitializer();
		initializer.initByName("estimate", false, "birthRate", birthRate, "speciesTree", transmissionTree, "geneTree",
				geneTree, "newick", initNewick);
		initializer.setID("SBI");


		FBDModel.initByName("origin", origin, "tree", transmissionTree, "birthRate", birthRate, "deathRate", deathRate,
				"samplingRate", samplingRate, "removalProbability", removalProbability);
		uniform3.initByName("upper", "1000.0");
		priorDist.initByName("x", origin, "distr", uniform3);
		
		List<Distribution> dist21 = new ArrayList<Distribution>();
		dist21.add(FBDModel);
		dist21.add(priorDist);
		prior.initByName("distribution", dist21);
		prior.setID("prior");

		// dist 2
		tnt.distribution.GeneTreeDistribution geneTreeDist = new tnt.distribution.GeneTreeDistribution();
		geneTreeDist.initByName("populationSizes", popSize, "origin", origin, "tau", tau,
				"birthRate", birthRate, "deathRate", deathRate, "samplingRate", samplingRate, "removalProbability",
				removalProbability,
				"geneTreeIntervals", intervals);
		List<Distribution> dist11 = new ArrayList<Distribution>();
		dist11.add(geneTreeDist);
		CompoundDistribution multiCoalescent = new CompoundDistribution();
		multiCoalescent.initByName("distribution", dist11);
		multiCoalescent.setID("multiCoalescent");

		// dist 3
		CompoundDistribution likelihood = new CompoundDistribution();

		TreeLikelihood treeLikelihood = new TreeLikelihood();
		
		SiteModel siteModel = new SiteModel();
		RealParameter mutationRate = new RealParameter("1.0");
		RealParameter shape = new RealParameter("1.0");
		RealParameter proportionInvariant = new RealParameter("0.0");
		JukesCantor substModel = new JukesCantor();
		mutationRate.initByName("estimate", false);
		shape.initByName("estimate", false);
		proportionInvariant.initByName("estimate", false, "lower", "0.0", "upper", "1.0");
		siteModel.initByName("mutationRate", mutationRate, "shape", shape, "proportionInvariant", proportionInvariant,
				"substModel", substModel);
		
		StrictClockModel branchRateModel = new StrictClockModel();
		RealParameter clockRate = new RealParameter("1.0");
		clockRate.initByName("estimate", false, "lower", "0.0");
		branchRateModel.initByName("clock.rate", clockRate);

		treeLikelihood.initByName("data", alignment, "tree", geneTree, "siteModel", siteModel, "branchRateModel",
				branchRateModel);

		List<Distribution> dist31 = new ArrayList<Distribution>();
		dist31.add(treeLikelihood);
		likelihood.initByName("distribution", dist31);
		likelihood.setID("likelihood");

		CompoundDistribution posterior = new CompoundDistribution();
		List<Distribution> distributions = new ArrayList<Distribution>();
		distributions.add(multiCoalescent);
		distributions.add(prior);
		distributions.add(likelihood);
		posterior.initByName("distribution", distributions);

		// Set up operators:
		// NodeReheight2
//		NodeReheight2 nodeReheight2 = new NodeReheight2();
//		nodeReheight2.initByName("tree", speciesTree, "geneTree", geneTreeDist, "taxonset", starTaxonSet, "weight",
//				"75.0");

		// CoordinatedUniform
		CoordinatedUniform coordinatedUniform = new CoordinatedUniform();
		coordinatedUniform.initByName("speciesTree", transmissionTree, "geneTree", geneTree, "weight", "30.0");

		// CoordinatedExponential
		CoordinatedExponential coordinatedExponential = new CoordinatedExponential();
		coordinatedExponential.initByName("speciesTree", transmissionTree, "geneTreeIntervals", intervals, "geneTree",
				geneTree, "weight", "30.0");

		// LeafToSampledAncestorJump
		LeafToSampledAncestorJump leafToSampleJump = new LeafToSampledAncestorJump();
		leafToSampleJump.initByName("tree", transmissionTree, "geneTreeIntervals", intervals,
				"weight", "75.0");

		// SAWilsonBalding
		SAWilsonBalding saWilsonBalding = new SAWilsonBalding();
		saWilsonBalding.initByName("tree", transmissionTree, "geneTreeIntervals", intervals,
				"weight", "50.0");

		// SAUniform
		SAUniform saUniform = new SAUniform();
		saUniform.initByName("tree", transmissionTree, "geneTree", geneTree, "geneTreeIntervals", intervals,
				"weight", "75.0");

//		// UpDownOperator_Species
//		UpDownOperator upDownOperatorSpecies = new UpDownOperator();
//		upDownOperatorSpecies.initByName("scaleFactor", "0.75", "weight", "6.0", "down", speciesTree, "down", geneTree);
		
//		// UpDownOperator_Species SA
//		UpDownOperator saUpDownOperatorSpecies = new UpDownOperator();
//		saUpDownOperatorSpecies.initByName("scaleFactor", "0.75", "weight", "6.0", "down", speciesTree);

//		// UpDownOperator_Gene
//		UpDownOperator upDownOperatorGene = new UpDownOperator();
//		upDownOperatorGene.initByName("scaleFactor", "0.95", "weight", "3.0","down", geneTree);

//		// ScaleOperator_GeneTree
//		ScaleOperator scaleOperatorGene = new ScaleOperator();
//		scaleOperatorGene.initByName("scaleFactor", "0.95", "weight", "3.0","tree", geneTree);
		
//		// ScaleOperator_GeneRoot
//		ScaleOperator scaleOperatorRoot = new ScaleOperator();
//		scaleOperatorRoot.initByName("rootOnly", true, "scaleFactor", "0.7", "weight", "3.0","tree", geneTree);

//		// Uniform
//		beast.evolution.operators.Uniform uniform = new beast.evolution.operators.Uniform();
//		uniform.initByName("weight", "15.0","tree", geneTree);

//		// SubtreeSlide
//		SubtreeSlide subtreeSlide = new SubtreeSlide();
//		subtreeSlide.initByName("size", "0.002", "tree", geneTree, "weight", "15.0");
//
//		// Exchange_Narrow
//		Exchange exchangeNarrow = new Exchange();
//		exchangeNarrow.initByName("tree", geneTree, "weight", "15.0");
//		
//		// Exchange_Wide
//		Exchange exchangeWide = new Exchange();
//		exchangeWide.initByName("isNarrow", false, "tree", geneTree, "weight", "15.0");

//		// WilsonBalding
//		WilsonBalding wilsonBalding = new WilsonBalding();
//		wilsonBalding.initByName("tree", geneTree, "weight", "15.0");

		// SPR
		SPROperatorChanged spr = new SPROperatorChanged();
		spr.initByName("tree", geneTree, "transmissionTree", transmissionTree, "geneTreeIntervals", intervals,
				"probBottleneck", "0.5", "rootAttachLambda", "0.5", "weight", "75.0");

		// Transmission_Attach
		TransmissionAttachChanged trAttach = new TransmissionAttachChanged();
		trAttach.initByName("tree", geneTree, "geneTreeIntervals", intervals,
				"rootAttachLambda", "0.5", "weight", "60.0");

		// Create_Mergers_Or_Reheight
		CreateMergersOrReheight cmor = new CreateMergersOrReheight();
		cmor.initByName("tree", geneTree, "mergerProb", "0.1", "transmissionTree", transmissionTree,
				"weight", "30.0");

		// Expand_Collapse
		ExpandCollapseOperator ec = new ExpandCollapseOperator();
		ec.initByName("tree", geneTree, "geneTreeIntervals", intervals,
				"weight", "30.0");

		// GeneExchange
		ExchangeOperator exchange = new ExchangeOperator();
		exchange.initByName("tree", geneTree, "geneTreeIntervals", intervals, "weight", "50.0");




//		// SAExchange_Wide
//		SAExchange saExchangeWide = new SAExchange();
//		saExchangeWide.initByName("isNarrow", false, "tree", speciesTree, "weight", "10.0");

//		// SAExchange_Narrow
//		SAExchange saExchangeNarrow = new SAExchange();
//		saExchangeNarrow.initByName("tree", speciesTree, "weight", "10.0");

//		// SAUniform
//		SAUniform saUniform = new SAUniform();
//		saUniform.initByName("tree", transmissionTree, "weight", "20.0");

//		// SAScaleOperator_Root
//		SAScaleOperator saScaleOperatorRoot = new SAScaleOperator();
//		saScaleOperatorRoot.initByName("rootOnly", true, "scaleFactor", "0.95", "tree", speciesTree, "weight", "1.0");
//
//		// SAScaleOperator_Root
//		SAScaleOperator saScaleOperatorTree = new SAScaleOperator();
//		saScaleOperatorTree.initByName("scaleFactor", "0.95", "tree", speciesTree, "weight", "3.0");

		Integer bunin = 1000;
		Integer chainLength = 50000000;
		Integer logEvery = 100;
		Integer statesLogged = (chainLength - bunin) / logEvery;

		// Set up logger:
		OrientedTreeLogger treeReport = new OrientedTreeLogger();
		treeReport.initByName("logEvery", logEvery.toString(),
				"burnin", bunin.toString(),
				"speciesTree", transmissionTree,
//				"geneTree", geneTreeDist,
				"log", transmissionTree,
				"silent", true);

		TransmissionTreeLogger tntLogger = new TransmissionTreeLogger();
		tntLogger.initByName("transmissionTree", transmissionTree);

		// Set up MCMC:
		MCMC mcmc = new MCMC();
		mcmc.initByName("chainLength", chainLength.toString(),
						"state", state, 
				"init", initializer,
						"distribution", posterior, 
//				"operator", nodeReheight2,
				"operator", coordinatedUniform,
				"operator", coordinatedExponential,
//				"operator", upDownOperatorSpecies,
//				"operator", upDownOperatorGene,
//				"operator", scaleOperatorGene,
//				"operator", scaleOperatorRoot,
//				"operator", uniform,
//				"operator", subtreeSlide,
//				"operator", exchangeNarrow,
//				"operator", exchangeWide,
//				"operator", originScaler,
				"operator", leafToSampleJump,
				"operator", saWilsonBalding,
//				"operator", saExchangeWide,
//				"operator", saExchangeNarrow,
				"operator", saUniform,
//				"operator", saScaleOperatorRoot,
//				"operator", saScaleOperatorTree,
				"operator", spr,
				"operator", trAttach,
				"operator", cmor,
				"operator", cmor,
				"operator", ec,
				"operator", exchange,
				"logger", treeReport);

		// Run MCMC:
		mcmc.run();

		int[][] frequencies = treeReport.getAnalysis();

		/*
		 * Eight Possible non-oriented topologies: ((3,2),1) ((3,1),2) (3,(2,1))
		 * ((3,2)1) (3,(2)1) ((3)1,2) ((3)2,1) (((3)2)1)
		 */

		// The 8 non-oriented topology frequencies have been calculated by Gavryushkina
		// et al.
		// (2014)
		double[] probs = new double[] { 0.778327, 0.043189, 0.043189, 0.078642, 0.006930, 0.006930, 0.038657,
				0.004135, };

		double tolerance = 0.025;
		double toleranceOriented = 0.025;
		double orientedFrequency = 0.25;

		for (int nonOrientedTopologyNr = 0; nonOrientedTopologyNr < 8; nonOrientedTopologyNr++) {
			int sumTopology = Arrays.stream(frequencies[nonOrientedTopologyNr]).sum();
			double probTopology = (double) sumTopology / (double) statesLogged;
			System.out.println("_____________________");
			System.out.println(probs[nonOrientedTopologyNr]);
			System.out.println("_____________________");
			System.out.println(probTopology);
			Assert.assertEquals(probs[nonOrientedTopologyNr], probTopology, tolerance);

			// For each non-oriented topology, there are four possible orientations
			// outputed.
			// We check for topologies with ancestral nodes separately below, because we do
			// not care for the orientation of Fake node children and therefore there are 2
			// possible orientations.
			// Last topology, where nodes 1 and 2 are direct ancestors can have only one
			// orientation and therefore is not validated.
			if (nonOrientedTopologyNr == 0 || nonOrientedTopologyNr == 1 || nonOrientedTopologyNr == 2) {
				for (int j = 0; j < 4; j++) {
					double frequency = (double) frequencies[nonOrientedTopologyNr][j] / (double) sumTopology;
					System.out.println(frequency);
					Assert.assertEquals(frequency, orientedFrequency, toleranceOriented);

				}
			}
			if (nonOrientedTopologyNr == 3) {
				for (int j = 0; j < 2; j++) {
					double frequency = (double) (frequencies[nonOrientedTopologyNr][j]
							+ frequencies[nonOrientedTopologyNr][j + 2]) / (double) sumTopology;
					System.out.println(frequency);
					Assert.assertEquals(frequency, 0.50, toleranceOriented);

				}
			}

			if (nonOrientedTopologyNr == 4 || nonOrientedTopologyNr == 5 || nonOrientedTopologyNr == 6) {
				for (int j = 0; j < 4; j += 2) {
					double frequency = (double) (frequencies[nonOrientedTopologyNr][j]
							+ frequencies[nonOrientedTopologyNr][j + 1]) / (double) sumTopology;
					System.out.println(frequency);
					Assert.assertEquals(frequency, 0.50, toleranceOriented);

				}
			}

		}
	}


	public class OrientedTreeLogger extends Logger {
		public Input<Integer> burninInput = new Input<Integer>("burnin",
				"Number of samples to skip (burn in)", Input.Validate.REQUIRED);

		public Input<Boolean> silentInput = new Input<Boolean>("silent",
				"Don't display final report.", false);

		final public Input<SpeciesTreeInterface> speciesTreeInput = new Input<>("speciesTree",
				"The species tree to be logged.", Validate.REQUIRED);

		final public Input<List<GeneTreeDistribution>> geneTreeInput = new Input<>("geneTree", "Gene tree within the species tree.",
				new ArrayList<>());

		private DecimalFormat df;


		SpeciesTreeInterface speciesTree;
		int m_nEvery = 1;
		int burnin;
		boolean silent = false;
		int[][] duplicate = new int[8][4];

		// There are 8 possible non-oriented topologies for a binary tree with two
		// ancestral samples and one extant sample. Each of those can be produced
		// in Beast2 in 4 different orientation.
		int[][] freq = new int[8][4];

		String number = "\\:(\\d+(\\.\\d+)?)(E\\-\\d+)?";

		// Beast2 will log all topologies in one of these three oriented groups.
		// The reason for this, is that a sampled ancestor node is added as a fake child
		// with edge of length 0;
		HashMap<Integer, Pattern> rx_1 = new HashMap<Integer, Pattern>() {
			{ // Group 1
				// (1,(2,3))
				put(1, Pattern.compile("\\(1" + number + "\\,\\(2" + number + "\\,3" + number + "\\)(.*)\\)(.*)"));
				// (1,(3,2))
				put(2, Pattern.compile("\\(1" + number + "\\,\\(3" + number + "\\,2" + number + "\\)(.*)\\)(.*)"));
				// ((2,3),1)
				put(3, Pattern.compile("\\(\\(2" + number + "\\,3" + number + "\\)(.*)\\,1" + number + "\\)(.*)"));
				// ((3,2),1)
				put(4, Pattern.compile("\\(\\(3" + number + "\\,2" + number + "\\)(.*)\\,1" + number + "\\)(.*)"));
			}
		};

		HashMap<Integer, Pattern> rx_2 = new HashMap<Integer, Pattern>() {
			{ // Group 2
				// (2, (1,3))
				put(1, Pattern.compile("\\(2" + number + "\\,\\(1" + number + "\\,3" + number + "\\)(.*)\\)(.*)"));
				// (2,(3,1))
				put(2, Pattern.compile("\\(2" + number + "\\,\\(3" + number + "\\,1" + number + "\\)(.*)\\)(.*)"));
				// ((1,3),2)
				put(3, Pattern.compile("\\(\\(1" + number + "\\,3" + number + "\\)(.*)\\,2" + number + "\\)(.*)"));
				// ((3,1),2)
				put(4, Pattern.compile("\\(\\(3" + number + "\\,1" + number + "\\)(.*)\\,2" + number + "\\)(.*)"));
			}
		};

		HashMap<Integer, Pattern> rx_3 = new HashMap<Integer, Pattern>() {
			{ // Group 3
				// (3,(1,2))
				put(1, Pattern.compile("\\(3" + number + "\\,\\(1" + number + "\\,2" + number + "\\)(.*)\\)(.*)"));
				// (3,(2,1))
				put(2, Pattern.compile("\\(3" + number + "\\,\\(2" + number + "\\,1" + number + "\\)(.*)\\)(.*)"));
				// ((2,1),3)
				put(3, Pattern.compile("\\(\\(2" + number + "\\,1" + number + "\\)(.*)\\,3" + number + "\\)(.*)"));
				// ((1,2),3)
				put(4, Pattern.compile("\\(\\(1" + number + "\\,2" + number + "\\)(.*)\\,3" + number + "\\)(.*)"));
			}
		};

		HashMap<Integer, HashMap<Integer, Pattern>> rx = new HashMap<Integer, HashMap<Integer, Pattern>>(){{
			put(1, rx_1);
			put(2, rx_2);
			put(3, rx_3);
		}};

		@Override
		public void initAndValidate() {

			List<BEASTObject> loggers = loggersInput.get();
			final int nLoggers = loggers.size();
			if (nLoggers == 0) {
				throw new IllegalArgumentException("Logger with nothing to log specified");
			}

			if (everyInput.get() != null)
				m_nEvery = everyInput.get();

			burnin = burninInput.get();

			if (silentInput.get() != null)
				silent = silentInput.get();

			speciesTree = speciesTreeInput.get();

		}

		@Override
		public void init() {
		}

		@Override
		public void log(long nSample) {

			if ((nSample % m_nEvery > 0) || nSample < burnin)
				return;

			SpeciesTreeInterface tree = (SpeciesTreeInterface) speciesTree.getCurrent();

			String newick = toNewick(tree.getRoot());
//			System.out.println(newick);

			boolean firstAncestor = false;
			boolean secondAncestor = false;
			Node first = null;
			List<Node> leaves = tree.getRoot().getAllLeafNodes();
			for (Node l : leaves) {
				if (l.getID().equals("1")) {
					first = l; 
					firstAncestor = l.isDirectAncestor() ? true : false;
				}
				if (l.getID().equals("2")) {
					secondAncestor = l.isDirectAncestor() ? true : false;
				}
			}
			

			for (int i = 0; i < freq.length; i++) {
				duplicate[i] = Arrays.copyOf(freq[i], freq[i].length);
			}


			if (!secondAncestor) {
				if (!firstAncestor) {
					// Non-oriented tree ((3,2),1) number 0, prob 77.8327
					// Non-oriented tree ((3,1),2) number 1, prob 4.3189
					// Non-oriented tree (3,(2,1)) number 2, prob 4.3189
					for (int nonOrientedTopologyNr=0; nonOrientedTopologyNr<3;nonOrientedTopologyNr++) {				
						for (int i = 1; i <= freq[nonOrientedTopologyNr].length; i++) {
							Matcher m = rx.get(nonOrientedTopologyNr + 1).get(i).matcher(newick);
							if (m.matches()) {
								freq[nonOrientedTopologyNr][i - 1] += 1;
							}
						}
					}
				} else if (firstAncestor) {
					// Non-oriented tree ((3,2)1) number 3:
					// node 1 is an ancestral node of both leaf nodes
					if (first.getParent().getNonDirectAncestorChild().getAllChildNodes().size() == 3) {
						int nonOrientedTopologyNr = 3;
						int orientedNr = 1;
						for (int i = 1; i <= freq[nonOrientedTopologyNr].length; i++) {
							Matcher m = rx.get(orientedNr).get(i).matcher(newick);
							if (m.matches()) {
								freq[nonOrientedTopologyNr][i - 1] += 1;
							}
						}
					}// Non-oriented tree (3,(2)1) number 4:
					// node 1 is an ancestral node of leaf node 2 
					else if (first.getParent().getNonDirectAncestorChild().getID() != null &&
							first.getParent().getNonDirectAncestorChild().getID().equals("2")) {
						int nonOrientedTopologyNr = 4;
						int orientedNr = 3;
						for (int i = 1; i <= freq[nonOrientedTopologyNr].length; i++) {
							if (i == 4) {
							}
							Matcher m = rx.get(orientedNr).get(i).matcher(newick);
							if (m.matches()) {
								freq[nonOrientedTopologyNr][i - 1] += 1;
							}
						}
					} // Non-oriented tree (2,(3)1) number 5:
					// node 1 is an ancestral node of leaf node 3
					else if (first.getParent().getNonDirectAncestorChild().getID() != null &&
							first.getParent().getNonDirectAncestorChild().getID().equals("3")) {
						int nonOrientedTopologyNr = 5;
						int orientedNr = 2;
						for (int i = 1; i <= freq[nonOrientedTopologyNr].length; i++) {
							Matcher m = rx.get(orientedNr).get(i).matcher(newick);
							if (m.matches()) {
								freq[nonOrientedTopologyNr][i - 1] += 1;
							}
						}
					}
				}
			}// Non-oriented tree (1,(3)2) number 6:
			// node 2 is an ancestral node of leaf node 3 
			else if (!firstAncestor) {
				int nonOrientedTopologyNr = 6;
				int orientedNr = 1;
				for (int i = 1; i <= freq[nonOrientedTopologyNr].length; i++) {
					Matcher m = rx.get(orientedNr).get(i).matcher(newick);
					if (m.matches()) {
						freq[nonOrientedTopologyNr][i - 1] += 1;
					}
				}
			} // Non-oriented tree (((3)2)1) number 7:
				// 1 is ancestral of 2 and 3, 2 is ancestral of 3
			else if (firstAncestor) {
				int nonOrientedTopologyNr = 7;
				int orientedNr = 1;
				for (int i = 1; i <= freq[nonOrientedTopologyNr].length; i++) {
					Matcher m = rx.get(orientedNr).get(i).matcher(newick);
					if (m.matches()) {
						freq[nonOrientedTopologyNr][i - 1] += 1;
					}
				}
			} else
				throw new IllegalArgumentException("Topology assignment fell through.");

			if (Arrays.equals(duplicate, freq)) {
				throw new IllegalArgumentException("Not Assigned: " + newick);
			}

			int sumBefore = 0;
			int sumAfter = 0;

			for (int nonOrientedTopologyNr = 0; nonOrientedTopologyNr < 8; nonOrientedTopologyNr++) {
				sumAfter += Arrays.stream(freq[nonOrientedTopologyNr]).sum();
				sumBefore += Arrays.stream(duplicate[nonOrientedTopologyNr]).sum();
			}

			if (sumAfter != sumBefore + 1) {
				throw new IllegalArgumentException("Tree assigned to more than one topology");
			}

		}

		@Override
		public void close() {
		}

		/**
		 * Obtain all frequencies
		 *
		 * @return topology frequencies.
		 */
		public int[][] getAnalysis() {
			return freq;
		}

		// Get Newick string, without sorting the nodes before
		String toNewick(Node node) {
			StringBuffer buf = new StringBuffer();
			if (node.getLeft() != null) {
				buf.append("(");
				buf.append(toNewick(node.getLeft()));
				if (node.getRight() != null) {
					buf.append(',');
					buf.append(toNewick(node.getRight()));
				}
				buf.append(")");
			}
			if (node.getID() == null) {
				buf.append(node.getNr() + 1);
			} else {
				buf.append(node.getID());
			}
			buf.append(":");

			double nodeLength;
			if (node.isRoot()) {
				nodeLength = getTreeHeight() - node.getHeight();
			} else {
				nodeLength = node.getLength();
			}
			appendDouble(buf, nodeLength);

			return buf.toString();
		}

		/**
		 * Appends a double to the given StringBuffer, formatting it using the private
		 * DecimalFormat instance, if the input 'dp' has been given a non-negative
		 * integer, otherwise just uses default formatting.
		 * 
		 * @param buf
		 * @param d
		 */
		private void appendDouble(StringBuffer buf, double d) {
			if (df == null) {
				buf.append(d);
			} else {
				buf.append(df.format(d));
			}
		}

		// uses the height of the tallest species or gene tree
		private double getTreeHeight() {
			double speciesTreeHeight = speciesTreeInput.get().getRoot().getHeight();

			for (GeneTreeDistribution gt : geneTreeInput.get()) {
				speciesTreeHeight = Double.max(speciesTreeHeight, gt.getRoot().getHeight());
			}

			return speciesTreeHeight;
		}

	}
	
}

