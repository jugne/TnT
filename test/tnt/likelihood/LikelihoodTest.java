package tnt.likelihood;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.math.util.MathUtils;
import org.junit.Test;

import beast.core.BEASTObject;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Logger;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.util.Randomizer;
import beast.util.TreeParser;
import feast.simulation.GPSimulator;
import tnt.simulator.SimulatedGeneTree;

public class LikelihoodTest {

	String speciesTreeNewick = "((t2_1:0.674025510498478,t7_1:0.674025510498478):1.02026915570423,((t6_1:0.441898780191674,(t11_1:0.321415277712137,t10_1:0.321415277712137):0.120483502479536):0.747685332950928,t22_1:0):0.504710553060102):2.13804639428371;";

	@Test
	public void topologyDistribution() throws Exception {

		Randomizer.setSeed(23255466);

		TaxonSet taxonSet = new TaxonSet();
		taxonSet.setID("taxonSet");
		Taxon t7_1 = new Taxon();
		t7_1.setID("t7_1");
		Taxon t10_1 = new Taxon();
		t10_1.setID("t10_1");
		Taxon t11_1 = new Taxon();
		t11_1.setID("t11_1");
		Taxon t6_1 = new Taxon();
		t6_1.setID("t6_1");
		Taxon t2_1 = new Taxon();
		t2_1.setID("t2_1");
		Taxon t22_1 = new Taxon();
		t22_1.setID("t22_1");

		taxonSet.initByName("taxon", t7_1, "taxon", t10_1, "taxon", t11_1, "taxon", t2_1, "taxon", t22_1, "taxon",
				t6_1);

		TreeParser transmissionTreeInput = new TreeParser();
		transmissionTreeInput.initByName("newick", speciesTreeNewick, "adjustTipHeights", false, "IsLabelledNewick",
				true);

		TraitSet sampleCounts = new TraitSet();
		sampleCounts.initByName("traitname", "sampleCounts", "taxa", taxonSet, "value", " t7_1=10,\n" +
				"t10_1=10,\n" +
				"t11_1=10,\n" +
				"t6_1=10,\n" +
				"t2_1=10,\n" +
				"t22_1=10");

		RealParameter populationSizes = new RealParameter("1.0");
		RealParameter bottleneckStrength = new RealParameter("4.0");

		RealParameter birthRate = new RealParameter("1.0");
		RealParameter deathRate = new RealParameter("0.4");
		RealParameter samplingRate = new RealParameter("0.3");

		SimulatedGeneTree geneTree = new SimulatedGeneTree();
		geneTree.setID("gene_tree_truth");
		geneTree.initByName("transmissionTreeInput", transmissionTreeInput,
				"sampleCounts", sampleCounts,
				"populationSizes", populationSizes,
				"bottleneckStrength", bottleneckStrength,
				"birthRate", birthRate,
				"deathRate", deathRate,
				"samplingRate", samplingRate);

		GeneTreeIntervals intervals = new GeneTreeIntervals();
		intervals.initByName("simulatedGeneTree", geneTree, "transmissionTreeInput", transmissionTreeInput);

		Integer runs = 100000;
		Integer burnin = 0;
		Integer logEvery = 1;

		LikelihoodDerivativeLogger log = new LikelihoodDerivativeLogger();
		log.initByName("logEvery", logEvery.toString(),
				"burnin", burnin.toString(),
				"geneTree", geneTree,
				"log", geneTree,
				"treeIntervals", intervals,
				"populationSizes", populationSizes,
				"bottleneckStrength", bottleneckStrength,
				"birthRate", birthRate,
				"deathRate", deathRate,
				"samplingRate", samplingRate);

		GPSimulator sim = new GPSimulator();
		sim.initByName("nSims", runs, "simulationObject", geneTree, "logger", log);

		sim.run();
		
		List<Double> derivatives = log.getAnalysis();
		double sum = derivatives.stream()
				.mapToDouble(a -> a)
				.sum();

		System.out.println("Here: " + sum / derivatives.size());


	}

	public class LikelihoodDerivativeLogger extends Logger {

		public Input<Integer> burninInput = new Input<Integer>("burnin",
				"Number of samples to skip (burn in)", Input.Validate.REQUIRED);

		final public Input<SimulatedGeneTree> geneTreeInput = new Input<>("geneTree",
				"Gene tree within the species tree.",
				Input.Validate.REQUIRED);

		final public Input<GeneTreeIntervals> treeIntervalsInput = new Input<>("treeIntervals",
				"Intervals for a phylogenetic beast tree",
				Input.Validate.REQUIRED);

		public Input<RealParameter> popSizesInput = new Input<RealParameter>("populationSizes",
				"Constant per-branch effective population sizes.", Validate.REQUIRED);


		public Input<RealParameter> bottleneckStrengthInput = new Input<RealParameter>("bottleneckStrength",
				"Strength of the bottleneck in scaled time", Validate.REQUIRED);


		public Input<RealParameter> birthRateInput = new Input<RealParameter>("birthRate", "Birth rate",
				Input.Validate.REQUIRED);

		public Input<Function> deathRateInput = new Input<Function>("deathRate", "Death rate", Input.Validate.REQUIRED);

		public Input<RealParameter> samplingRateInput = new Input<RealParameter>("samplingRate",
				"Sampling rate per individual", Input.Validate.REQUIRED);



		SimulatedGeneTree geneTree;

		GeneTreeIntervals intervals;
		int m_nEvery = 1;
		int burnin;
		double Ne;
		double tau;
		double lambda;
		double mu;
		double psi;

		int s;

		List<Double> derivatives = new ArrayList<Double>();
//		int[] endLineages = new int[6];

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
			geneTree = geneTreeInput.get();
			intervals = treeIntervalsInput.get();

			Ne = popSizesInput.get().getArrayValue(0);
			tau = bottleneckStrengthInput.get().getArrayValue(0);
			lambda = birthRateInput.get().getArrayValue(0);
			mu = deathRateInput.get().getArrayValue(0);
			psi = samplingRateInput.get().getArrayValue(0);

			s = 0;

		}

		@Override
		public void init() {
		}

		@Override
		public void close() {
		}

		@Override
		public void log(long nSample) {

			if ((nSample % m_nEvery > 0) || nSample < burnin)
				return;

			SimulatedGeneTree tree = (SimulatedGeneTree) geneTree.getCurrent();
			intervals.requiresRecalculation();
			HashMap<Node, List<GeneTreeEvent>> eventList = intervals.getGeneTreeEventList();
			double h = 0.000001;
			double derivative = 0;

//			System.out.println(geneTree.getRoot().toNewick());

			List<Node> trNodeList = new ArrayList<Node>(eventList.keySet());
			trNodeList = trNodeList.stream().sorted(Comparator.comparingDouble(n -> n.getHeight()))
					.collect(Collectors.toList());


			Node trNode = trNodeList.get(5);
			boolean donor = (trNode.getParent().getChild(0) == trNode && !trNode.getParent().isFake());

			List<GeneTreeEvent> localEventList = eventList.get(trNode).stream()
					.sorted(Comparator.comparingDouble(e -> e.time))
					.collect(Collectors.toList());
			GeneTreeEvent prevEvent = new GeneTreeEvent();
			prevEvent.time = trNode.getHeight();
			for (GeneTreeEvent event : localEventList) {
//				if (trNode.getParent().getHeight() - 0.001 <= event.time)
//					break;

				// Contribution from every interval, except the last
				if (prevEvent.time < event.time) {
					derivative += logInterval(event, prevEvent, Ne + h);
					derivative -= logInterval(event, prevEvent, Ne);
				}

				switch (event.type) {
				case SAMPLE:
					break;
				case BIFURCATION:
					// If event is recorded at time of transmission and
					// transmission tree branch is recipient,
					// log as transmission event.
					// This means no multiplication with rate \lambda*P_0
					if (!trNode.isRoot() && !donor && event.time == trNode.getParent().getHeight()) {
						derivative += logTrans(event, prevEvent, Ne + h);
						derivative -= logTrans(event, prevEvent, Ne);
						break;
					}
					derivative += logBif(event, prevEvent, Ne + h);
					derivative -= logBif(event, prevEvent, Ne);
					break;
				case MULTIFURCATION:
					// To see how many multiple mergers we have so far
					if (event.multiCoalSize.size() > 1)
						System.out.println(++s);

					// If event is recorded at time of transmission and
					// transmission tree branch is recipient,
					// log as transmission event.
					// This means no multiplication with rate \lambda*P_0
					if (!trNode.isRoot() && !donor && event.time == trNode.getParent().getHeight()) {
						derivative += logTrans(event, prevEvent, Ne + h);
						derivative -= logTrans(event, prevEvent, Ne);
							break;
						}
					derivative += logMultif(event, prevEvent, Ne + h);
					derivative -= logMultif(event, prevEvent, Ne);
				}
				prevEvent = event;
			}

			if (!trNode.isRoot() && prevEvent.time < trNode.getParent().getHeight()) {
						GeneTreeEvent mockEvent = new GeneTreeEvent();
				mockEvent.time = trNode.getParent().getHeight();// - 0.001;
				mockEvent.lineages = prevEvent.lineages;
				derivative += logInterval(mockEvent, prevEvent, Ne + h);
				derivative -= logInterval(mockEvent, prevEvent, Ne);
//				if (!donor) {
//					derivative += logTrans(mockEvent, prevEvent, Ne + h);
//					derivative -= logTrans(mockEvent, prevEvent, Ne);
//				}
			}
			derivatives.add(derivative / h);
		}

		public List<Double> getAnalysis() {
//			double[] normal = new double[endLineages.length];
//			for (int a = 0; a < endLineages.length; a++) {
//				normal[a] = endLineages[a] / (double) IntStream.of(endLineages).sum();
//			}
//			System.out.println("Frequencies: " + Arrays.toString(normal));

			double derivatives_100 = derivatives.subList(0, 100).stream()
					.mapToDouble(a -> a)
					.sum();
			System.out.println("100 runs: " + derivatives_100 / 100);
			double derivatives_1000 = derivatives.subList(0, 1000).stream()
					.mapToDouble(a -> a)
					.sum();
			System.out.println("1000 runs: " + derivatives_1000 / 1000);
			double derivatives_10000 = derivatives.subList(0, 10000).stream()
					.mapToDouble(a -> a)
					.sum();
			System.out.println("1000 runs: " + derivatives_10000 / 10000);
			double derivatives_100000 = derivatives.subList(0, 100000).stream()
					.mapToDouble(a -> a)
					.sum();
			System.out.println("10000 runs: " + derivatives_100000 / 100000);

			return derivatives;
		}


		private double logBif(GeneTreeEvent event, GeneTreeEvent prevEvent, double Ne) {
			double ans = 0.0;
			ans += 1.0 / Ne;
			ans += (1.0 / waysToCoal(prevEvent.lineages, event.lineages))
					* gUp(prevEvent.lineages, event.lineages, tau, Ne)
					* lambda * P_0(event.time);

			return Math.log(ans);
		}

		private double logMultif(GeneTreeEvent event, GeneTreeEvent prevEvent, double Ne) {
			double ans = 0.0;
			double mult = 1.0;
			int sum = event.multiCoalSize.stream().mapToInt(Integer::intValue)
					.sum();
			int n_histories = event.multiCoalSize.size();
			for (int s=0; s<n_histories; s++) {
				mult *= waysToCoal(event.multiCoalSize.get(s), 1);
				// W factor from NOAH A. ROSENBERG
				mult *= binomialInt(sum - n_histories, event.multiCoalSize.get(s) - 1);
				sum -= (event.multiCoalSize.get(s) - 1);
			}
			ans += (1.0 / waysToCoal(prevEvent.lineages, event.lineages)) * mult
					* gUp(prevEvent.lineages, event.lineages, tau, Ne)
					* lambda * P_0(event.time);

			return Math.log(ans);
		}

		private double logTrans(GeneTreeEvent event, GeneTreeEvent prevEvent, double Ne) {
			double ans = 0.0;
			double mult = 1.0;
			int sum = event.multiCoalSize.stream().mapToInt(Integer::intValue)
					.sum();
			int n_histories = event.multiCoalSize.size();
			for (int s = 0; s < n_histories; s++) {
				mult *= waysToCoal(event.multiCoalSize.get(s), 1);
				// W factor from NOAH A. ROSENBERG
				mult *= binomialInt(sum - n_histories, event.multiCoalSize.get(s) - 1);
				sum -= (event.multiCoalSize.get(s) - 1);
			}
			ans += (1.0 / waysToCoal(prevEvent.lineages, event.lineages)) * mult
					* gUp(prevEvent.lineages, event.lineages, tau, Ne);

			return Math.log(ans);
		}

		private double logInterval(GeneTreeEvent event, GeneTreeEvent prevEvent, double Ne) {
			double ans = 0.0;
			ans -= prevEvent.lineages * (prevEvent.lineages - 1) * 0.5 * (1.0 / Ne) * (event.time - prevEvent.time);

			double sum = 0;
			for (int j = 1; j < prevEvent.lineages; j++) {
				sum += gUp(prevEvent.lineages, j, tau, Ne);

			}

			ans -= sum * lambda
					* integralP_0(prevEvent.time, event.time, lambda, mu, psi);

			return ans;
		}

		private double gUp(int i, int j, double tau, double Ne) {
			double ans = 0.0;
			for (int k = j; k <= i; k++) {
				ans += (2 * k - 1) * Math.pow(-1, k - j) * f_1(j, k - 1) * f_2(i, k)
						* Math.exp(-(k * (k - 1) * tau * 0.5) / Ne) /
						(MathUtils.factorial(j) * MathUtils.factorial(k - j) * f_1(i, k));
			}
			
			return ans;
		}

		private double P_0(double t) {

			double c_1 = Math.abs(Math.sqrt(Math.pow((lambda - mu - psi), 2) + 4 * lambda * psi));
			double c_2 = -(lambda - mu - 2 * lambda - psi) / c_1;

			
			return (lambda + mu + psi
					+ c_1 * ((Math.exp(-c_1 * t) * (1 - c_2) - (1 + c_2))
							/ (Math.exp(-c_1 * t) * (1 - c_2) + (1 + c_2))))
					/ (2.0 * lambda);
		}

		private double integralP_0(double t_0, double t_1, double lambda, double mu, double psi) {
			double c_1 = Math.abs(Math.sqrt(Math.pow((lambda - mu - psi), 2) + 4 * lambda * psi));
			double c_2 = -(lambda - mu - 2 * lambda - psi) / c_1;

			double ans = (1.0 / (2.0 * lambda)) * ((t_1 - t_0) * (mu + psi + lambda - c_1) + 2.0 * Math
					.log(((c_2 - 1) * Math.exp(-c_1 * t_0) - c_2 - 1) / ((c_2 - 1) * Math.exp(-c_1 * t_1) - c_2 - 1)));

			return ans;
		}

		private double waysToCoal(int i, int j) {
			double ans = MathUtils.factorial(i);
			ans *= MathUtils.factorial(i - 1);
			ans /= Math.pow(2, i - j);
			ans /= MathUtils.factorial(j);
			ans /= MathUtils.factorial(j - 1);

			return ans;
		}

		private double f_1(int x, int y) {
			double ans = 1.0;
			for (int i = 0; i < y; i++) {
				ans *= (x + i);
			}

			return ans;
		}

		private double f_2(int x, int y) {
			double ans = 1.0;
			for (int i = 0; i < y; i++) {
				ans *= (x - i);
			}

			return ans;
		}

		private long binomialInt(int n, int k) {
			if (k > n - k)
				k = n - k;

			long binom = 1;
			for (int i = 1; i <= k; i++)
				binom = binom * (n + 1 - i) / i;
			return binom;
		}

	}

}
