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
		sampleCounts.initByName("traitname", "sampleCounts", "taxa", taxonSet, "value", " t7_1=4,\n" +
				"t10_1=1,\n" +
				"t11_1=1,\n" +
				"t6_1=1,\n" +
				"t2_1=1,\n" +
				"t22_1=1");

		RealParameter populationSizes = new RealParameter("1.0");
		RealParameter bottleneckStrength = new RealParameter("1.0");

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

		Integer runs = 10000000;
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
		int kl;

		List<Double> derivatives = new ArrayList<Double>();
		List<Double> derivatives_2 = new ArrayList<Double>();
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
			kl = 1;

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
			double h = 0.0001;
			double derivative = 0;
			double derivative_1 = 1;
			double derivative_2 = 1;

//			System.out.println(geneTree.getRoot().toNewick());

			List<Node> trNodeList = new ArrayList<Node>(eventList.keySet());
			trNodeList = trNodeList.stream().sorted(Comparator.comparingDouble(n -> n.getHeight()))
					.collect(Collectors.toList());


			Node trNode = trNodeList.get(0);
			boolean donor = (trNode.getParent().getChild(0) == trNode && !trNode.getParent().isFake());

			List<GeneTreeEvent> localEventList = eventList.get(trNode).stream()
					.sorted(Comparator.comparingDouble(e -> e.time))
					.collect(Collectors.toList());
			GeneTreeEvent prevEvent = new GeneTreeEvent();
			prevEvent.time = trNode.getHeight();
			for (GeneTreeEvent event : localEventList) {
				if (event.lineages == 0)
					System.out.println(geneTree.getRoot().toNewick());
//				if (trNode.getParent().getHeight() - 0.001 <= event.time)
//					break;

				// Contribution from every interval, except the last
				if (prevEvent.time < event.time) {
//					derivative_1 *= logInterval(event, prevEvent, Ne + h);
//					derivative_2 *= logInterval(event, prevEvent, Ne - h);
					derivative += derivativeNeLogInt(event, prevEvent);
//					derivative += Math.log(logInterval(event, prevEvent, Ne + h));
//					derivative -= Math.log(logInterval(event, prevEvent, Ne - h));
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
//						derivative_1 *= logTrans(event, prevEvent, Ne + h);
//						derivative_2 *= logTrans(event, prevEvent, Ne - h);
						derivative += derivativeNeLogMulti(event, prevEvent);
//						derivative += Math.log(logTrans(event, prevEvent, Ne + h));
//						derivative -= Math.log(logTrans(event, prevEvent, Ne - h));
						break;
					}
//					derivative_1 *= logBif(event, prevEvent, Ne + h);
//					derivative_2 *= logBif(event, prevEvent, Ne - h);
					derivative += derivativeNeLogBif(event, prevEvent);
//					derivative += Math.log(logBif(event, prevEvent, Ne + h));
//					derivative -= Math.log(logBif(event, prevEvent, Ne - h));
					break;
				case MULTIFURCATION:
					// To see how many multiple mergers we have so far
//					if (event.multiCoalSize.size() > 1) {
//						System.out.println(++s);
//						if (s == 3)
//							System.out.println();
////						if (s == 6373)
////							System.out.println(geneTree.getRoot().toNewick());
//					}

					// If event is recorded at time of transmission and
					// transmission tree branch is recipient,
					// log as transmission event.
					// This means no multiplication with rate \lambda*P_0
					if (!trNode.isRoot() && !donor && event.time == trNode.getParent().getHeight()) {
//						derivative_1 *= logTrans(event, prevEvent, Ne + h);
//						derivative_2 *= logTrans(event, prevEvent, Ne - h);
						derivative += derivativeNeLogMulti(event, prevEvent);
//						derivative += Math.log(logTrans(event, prevEvent, Ne + h));
//						derivative -= Math.log(logTrans(event, prevEvent, Ne - h));
							break;
						}
//					derivative_1 *= logMultif(event, prevEvent, Ne + h);
//					derivative_2 *= logMultif(event, prevEvent, Ne - h);
					derivative += derivativeNeLogMulti(event, prevEvent);
//					derivative += Math.log(logMultif(event, prevEvent, Ne + h));
//					derivative -= Math.log(logMultif(event, prevEvent, Ne - h));
				}
				prevEvent = event;
			}

			if (!trNode.isRoot() && prevEvent.time < trNode.getParent().getHeight()) {
						GeneTreeEvent mockEvent = new GeneTreeEvent();
				mockEvent.time = trNode.getParent().getHeight();// - 0.001;
				mockEvent.lineages = prevEvent.lineages;
//				derivative_1 *= logInterval(mockEvent, prevEvent, Ne + h);
//				derivative_2 *= logInterval(mockEvent, prevEvent, Ne - h);
				derivative += derivativeNeLogInt(mockEvent, prevEvent);
//				derivative += Math.log(logInterval(mockEvent, prevEvent, Ne + h));
//				derivative -= Math.log(logInterval(mockEvent, prevEvent, Ne - h));
				if (!donor) {
//					derivative_1 *= logTrans(mockEvent, prevEvent, Ne + h);
//					derivative_2 *= logTrans(mockEvent, prevEvent, Ne - h);
					derivative += derivativeNeLogMulti(mockEvent, prevEvent);
//					derivative += Math.log(logTrans(mockEvent, prevEvent, Ne + h));
//					derivative -= Math.log(logTrans(mockEvent, prevEvent, Ne - h));
				}
			}
//			double m = Math.log(derivative_1) - Math.log(derivative_2);
//			m /= (2 * h);
////			derivatives.add(derivative / (2 * h));
//			derivatives.add(m);

			derivatives_2.add(derivative);

			if (derivatives_2.size() == 100) {
//				double derivatives_100 = derivatives.subList(0, 100).stream()
//						.mapToDouble(a -> a)
//						.sum();

				double derivatives_100_2 = derivatives_2.subList(0, 100).stream()
						.mapToDouble(a -> a)
						.sum();

//				System.out.println("100 runs: " + derivatives_100 / 100);
				System.out.println("100 runs: " + derivatives_100_2 / 100);
			}

			if (derivatives_2.size() == 1000) {
//				double derivatives_1000 = derivatives.subList(0, 1000).stream()
//						.mapToDouble(a -> a)
//						.sum();

				double derivatives_1000_2 = derivatives_2.subList(0, 1000).stream()
						.mapToDouble(a -> a)
						.sum();

//				System.out.println("1000 runs: " + derivatives_1000 / 1000);
				System.out.println("1000 runs: " + derivatives_1000_2 / 1000);
			}

			if (derivatives_2.size() / 10000.0 == derivatives_2.size() / 10000 && derivatives_2.size() / 10000 >= 1
					&& derivatives_2.size() / 10000 < 10) {
//				double derivatives_10000 = derivatives.subList(0, kl * 10000).stream()
//						.mapToDouble(a -> a)
//						.sum();

				double derivatives_10000_2 = derivatives_2.subList(0, kl * 10000).stream()
						.mapToDouble(a -> a)
						.sum();

//				System.out.println("10000 runs: " + derivatives_10000 / (kl * 10000));
				System.out.println(kl + "x 10000 runs: " + derivatives_10000_2 / (kl * 10000));
				kl++;
			}

			if (derivatives_2.size() == 100000) {
//				double derivatives_100000 = derivatives.subList(0, 100000).stream()
//						.mapToDouble(a -> a)
//						.sum();
				double derivatives_100000_2 = derivatives_2.subList(0, 100000).stream()
						.mapToDouble(a -> a)
						.sum();

//				System.out.println("100000 runs: " + derivatives_100000 / 100000);
				System.out.println("100000 runs: " + derivatives_100000_2 / 100000);
			}

			if (derivatives_2.size() == 1000000) {
//				double derivatives_1000000 = derivatives.subList(0, 1000000).stream()
//						.mapToDouble(a -> a)
//						.sum();
				double derivatives_1000000_2 = derivatives_2.subList(0, 1000000).stream()
						.mapToDouble(a -> a)
						.sum();

//				System.out.println("1000000 runs: " + derivatives_1000000 / 1000000);
				System.out.println("1000000 runs: " + derivatives_1000000_2 / 1000000);
			}

			if (derivatives_2.size() == 10000000) {
//				double derivatives_10000000 = derivatives.subList(0, 1000000).stream()
//						.mapToDouble(a -> a)
//						.sum();
				double derivatives_10000000_2 = derivatives_2.subList(0, 1000000).stream()
						.mapToDouble(a -> a)
						.sum();

//				System.out.println("10000000 runs: " + derivatives_10000000 / 10000000);
				System.out.println("10000000 runs: " + derivatives_10000000_2 / 10000000);
			}

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

			double derivatives_100_2 = derivatives_2.subList(0, 100).stream()
					.mapToDouble(a -> a)
					.sum();

			System.out.println("100 runs: " + derivatives_100 / 100);
			System.out.println("100_2 runs: " + derivatives_100_2 / 100);

			double derivatives_1000 = derivatives.subList(0, 1000).stream()
					.mapToDouble(a -> a)
					.sum();

			double derivatives_1000_2 = derivatives_2.subList(0, 1000).stream()
					.mapToDouble(a -> a)
					.sum();

			System.out.println("1000 runs: " + derivatives_1000 / 1000);
			System.out.println("1000_2 runs: " + derivatives_1000_2 / 1000);

			double derivatives_10000 = derivatives.subList(0, 10000).stream()
					.mapToDouble(a -> a)
					.sum();

			double derivatives_10000_2 = derivatives_2.subList(0, 10000).stream()
					.mapToDouble(a -> a)
					.sum();

			System.out.println("10000 runs: " + derivatives_10000 / 10000);
			System.out.println("10000_2 runs: " + derivatives_10000_2 / 10000);

			double derivatives_100000 = derivatives.subList(0, 100000).stream()
					.mapToDouble(a -> a)
					.sum();
			double derivatives_100000_2 = derivatives_2.subList(0, 100000).stream()
					.mapToDouble(a -> a)
					.sum();

			System.out.println("100000 runs: " + derivatives_100000 / 100000);
			System.out.println("100000_2 runs: " + derivatives_100000_2 / 100000);

			double derivatives_1000000 = derivatives.subList(0, 1000000).stream()
					.mapToDouble(a -> a)
					.sum();
			System.out.println("1000000 runs: " + derivatives_1000000 / 1000000);

			return derivatives;
		}


		private double logBif(GeneTreeEvent event, GeneTreeEvent prevEvent, double Ne) {
			double ans = 0.0;
			ans += 1.0 / Ne;
			ans += (1.0 / waysToCoal(prevEvent.lineages, event.lineages))
					* gUp(prevEvent.lineages, event.lineages, tau, Ne)
					* lambda * P_0(event.time);


			return ans;
//			return Math.log(ans);
		}

		private double logMultif(GeneTreeEvent event, GeneTreeEvent prevEvent, double Ne) {
			double ans;
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
			ans = (1.0 / waysToCoal(prevEvent.lineages, event.lineages)) * mult
					* gUp(prevEvent.lineages, event.lineages, tau, Ne)
					* lambda * P_0(event.time);

			return ans;
//			return Math.log(ans);
		}

		private double logTrans(GeneTreeEvent event, GeneTreeEvent prevEvent, double Ne) {
			double ans;
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
			ans = (1.0 / waysToCoal(prevEvent.lineages, event.lineages)) * mult
					* gUp(prevEvent.lineages, event.lineages, tau, Ne);

			return ans;
//			return Math.log(ans);
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

			return Math.exp(ans);
//			return ans;
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

		private double derivativeNeLogMulti(GeneTreeEvent event, GeneTreeEvent prevEvent) {
			int i = prevEvent.lineages;
			int j = event.lineages;
			double up = 0.0;
			double down = 0.0;
			for (int k = j; k <= i; k++) {
				double coef = (((2 * k - 1) * Math.pow(-1, k - j) * f_1(j, k - 1) * f_2(i, k)) /
						(MathUtils.factorial(j) * MathUtils.factorial(k - j) * f_1(i, k)));
				up += coef * (k * (k - 1) / (2 * Math.pow(Ne, 2))) * Math.exp(-(k * (k - 1) * tau * 0.5) / Ne);
				down += coef * Math.exp(-(k * (k - 1) * tau * 0.5) / Ne);
			}

			return up / down;
		}

		private double derivativeNeLogBif(GeneTreeEvent event, GeneTreeEvent prevEvent) {
			int i = prevEvent.lineages;
			int j = event.lineages;
			double up = 0.0;
			double down = 0.0;
			for (int k = j; k <= i; k++) {
				double coef = (((2 * k - 1) * Math.pow(-1, k - j) * f_1(j, k - 1) * f_2(i, k)) /
						(MathUtils.factorial(j) * MathUtils.factorial(k - j) * f_1(i, k)));
				up += coef * (k * (k - 1) / (2 * Math.pow(Ne, 2))) * Math.exp(-(k * (k - 1) * tau * 0.5) / Ne);
				down += coef * Math.exp(-(k * (k - 1) * tau * 0.5) / Ne);
			}
			

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
			
			up = ((-1 / Math.pow(Ne, 2))
					+ ((1.0 / waysToCoal(prevEvent.lineages, event.lineages)) * mult * lambda * P_0(event.time) *
							up));

			down = ((1 / Ne) +
					((1.0 / waysToCoal(prevEvent.lineages, event.lineages)) * mult * lambda * P_0(event.time) *
							down));

			return up / down;
		}

		private double derivativeNeLogInt(GeneTreeEvent event, GeneTreeEvent prevEvent) {
			int l = prevEvent.lineages;
			double t_1 = event.time;
			double t_0 = prevEvent.time;
			
			double ans = ((l * (l - 1)) * (t_1 - t_0) / (2 * Math.pow(Ne, 2)))
					+ (lambda * integralP_0(prevEvent.time, event.time, lambda, mu, psi)
							* ((2 * l - 1) * f_1(l, l - 1) * f_2(l, l) * l * (l - 1) * tau
									* Math.exp(-(l * (l - 1) * tau * 0.5) / Ne))
							/ (2 * MathUtils.factorial(l) * f_1(l, l) * Math.pow(Ne, 2)));
			
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
