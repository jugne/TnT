package tnt.likelihood;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.math.util.MathUtils;

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

// Copies intro and fome helper functions from src/coalre/networkannotator/ReassortmentAnnotator.java 
// in CoalRe package by N.F. Müller
public class LikelihoodTest {

	static String speciesTreeNewick = "((t2_1:0.674025510498478,t7_1:0.674025510498478):1.02026915570423,((t6_1:0.441898780191674,(t11_1:0.321415277712137,t10_1:0.321415277712137):0.120483502479536):0.747685332950928,t22_1:0):0.504710553060102):2.13804639428371;";

	public static String helpMessage = "LikelihoodTest - iteratively samples from gene tree distribution, given the fixed transmission tree.\n"
			+ "Then tests, if Expected Score statistic is zero (assimptotically).\n"
			+ " Test is done for one recipient bramnch of transmission tree.\n"
			+ "\n"
			+ "Option                   Description\n"
			+ "--------------------------------------------------------------\n"
			+ "-runs                    Choose how many iterations to run. (Default 1E7)\n"
			+ "-calcEvery				Choose how often to calculate score statistic. Cannot be larger than runs. (Default 1E2)\n"
			+ "-Ne						Effective population size. (Default 1.0)\n"
			+ "-tau						Bottleneck duration/streangth. (Default 1.0)\n"
			+ "-seed					If set, the test will be reproducible excatly with the same seed.\n"
			+ "-nGeneSamples			NUmber of samples per transmission tree branch.\n"
			+ "\n"
			+ "Output is written to a file named 'likelihood_test.txt'.";

//	@Test
//	public void topologyDistribution() throws Exception {
	public static void main(String[] args) throws Exception {

		Integer runs = 10000000;
		Integer burnin = 0;
		Integer logEvery = 1;
		Integer calcEvery = 10;
		Long seed = Randomizer.nextLong();

		Double Ne = 1.0;
		Double tau = 1.0;
		Integer nGeneSamples = 2;

		int i = 0;
		if (args.length > 0) {
			while (args[i].startsWith("-")) {
				switch (args[i]) {
				case "-help":
					printUsageAndExit();
					break;
				case "-runs":
					if (args.length <= i + 1)
						printUsageAndError("-runs must be followed by a number");

					try {
						runs = Integer.parseInt(args[i + 1]);
					} catch (NumberFormatException e) {
						printUsageAndError("Error parsing runs number.");
					}

					if (runs <= 0) {
						printUsageAndError("Runs must be > 0.");
					}

					i += 1;
					break;
				case "-calcEvery":
					if (args.length <= i + 1)
						printUsageAndError("-calcEvery must be followed by a number");

					try {
						calcEvery = Integer.parseInt(args[i + 1]);
					} catch (NumberFormatException e) {
						printUsageAndError("Error parsing calcEvery.");
					}

					if (calcEvery <= 0) {
						printUsageAndError("CalcEvery must be > 0.");
					}

					i += 1;
					break;
				case "-Ne":
					if (args.length <= i + 1)
						printUsageAndError("-Ne must be followed by a number");

					try {
						Ne = Double.parseDouble(args[i + 1]);
					} catch (NumberFormatException e) {
						printUsageAndError("Error parsing Ne number.");
					}

					if (Ne <= 0) {
						printUsageAndError("Ne must be > 0.");
					}

					i += 1;
					break;
				case "-tau":
					if (args.length <= i + 1)
						printUsageAndError("-tau must be followed by a number");

					try {
						tau = Double.parseDouble(args[i + 1]);
					} catch (NumberFormatException e) {
						printUsageAndError("Error parsing bottleneck strength.");
					}

					if (tau <= 0) {
						printUsageAndError("Bottleneck strength must be > 0.");
					}

					i += 1;
					break;

				case "-seed":
					if (args.length <= i + 1)
						printUsageAndError("-seed must be followed by a number");

					try {
						seed = Long.parseLong(args[i + 1]);
					} catch (NumberFormatException e) {
						printUsageAndError("Error parsing seed number.");
					}

					if (seed <= 0) {
						printUsageAndError("seed must be > 0.");
					}

					i += 1;
					break;
				case "-nGeneSamples":
					if (args.length <= i + 1)
						printUsageAndError("-seed must be followed by a number");

					try {
						nGeneSamples = Integer.parseInt(args[i + 1]);
					} catch (NumberFormatException e) {
						printUsageAndError("Error parsing nGeneSamples number.");
					}

					if (nGeneSamples <= 0) {
						printUsageAndError("nGeneSamples must be > 0.");
					}

					i += 1;
					break;

				default:
					printUsageAndError("Unrecognised command line option '" + args[i] + "'.");
				}

				i += 1;
			}
		}


		Randomizer.setSeed(seed);

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
		sampleCounts.initByName("traitname", "sampleCounts", "taxa", taxonSet, "value",
				" t7_1=" + nGeneSamples + ",\n" +
				"t10_1=1,\n" +
				"t11_1=1,\n" +
				"t6_1=1,\n" +
				"t2_1=1,\n" +
				"t22_1=1");

		RealParameter populationSizes = new RealParameter(Ne.toString());
		RealParameter bottleneckStrength = new RealParameter(tau.toString());

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
				"samplingRate", samplingRate,
				"calcEvery", calcEvery);

		GPSimulator sim = new GPSimulator();
		sim.initByName("nSims", runs, "simulationObject", geneTree, "logger", log);
		sim.run();
		log.getAnalysis();

		System.out.println("Done");


	}

	public static class LikelihoodDerivativeLogger extends Logger {

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

		public Input<Integer> calcEveryInput = new Input<Integer>("calcEvery",
				"How often to calculate score and record empirical expectation", Input.Validate.REQUIRED);


		SimulatedGeneTree geneTree;
		PrintStream ps;
		GeneTreeIntervals intervals;
		int m_nEvery = 1;
		int calcEvery;
		int burnin;
		double Ne;
		double tau;
		double lambda;
		double mu;
		double psi;

		int kl;

		double sumDerivativesApprox;
		double sumDerivativesExact;

		List<Double> derivativesApprox = new ArrayList<Double>();
		List<Double> derivativesExact = new ArrayList<Double>();

		@Override
		public void initAndValidate() {

			List<BEASTObject> loggers = loggersInput.get();
			final int nLoggers = loggers.size();
			if (nLoggers == 0) {
				throw new IllegalArgumentException("Logger with nothing to log specified");
			}

			if (everyInput.get() != null)
				m_nEvery = everyInput.get();

			calcEvery = calcEveryInput.get();
			burnin = burninInput.get();
			geneTree = geneTreeInput.get();
			intervals = treeIntervalsInput.get();

			Ne = popSizesInput.get().getArrayValue(0);
			tau = bottleneckStrengthInput.get().getArrayValue(0);
			lambda = birthRateInput.get().getArrayValue(0);
			mu = deathRateInput.get().getArrayValue(0);
			psi = samplingRateInput.get().getArrayValue(0);

			kl = 1;
			
			sumDerivativesApprox = 0;
			sumDerivativesExact = 0;

			try {
				ps = new PrintStream("likelihood_test.txt");
				ps.print("n_run" + "\t" + "expected_score_approx" + "\t" + "expected_score_exact");
				ps.print("\n");
			} catch (FileNotFoundException e) {
				e.printStackTrace();
				System.exit(0);
			}

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

			double derivative = 0;

			// delta for approximate derivative
			double h = 0.0001;
			double derivativeApprox_1 = 1;
			double derivativeApprox_2 = 1;

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

				// Check that there are no events that leave non-positive number of lineages
				if (event.lineages == 0) {
					System.out.println("Zero lineages after event!!");
					System.exit(1);
				}

				// Contribution from every interval, except the last
				if (prevEvent.time < event.time) {
					derivativeApprox_1 *= interval(event, prevEvent, Ne + h);
					derivativeApprox_2 *= interval(event, prevEvent, Ne - h);
					derivative += derivativeNeLogInt(event, prevEvent);
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
						derivativeApprox_1 *= knownTransmission(event, prevEvent, Ne + h);
						derivativeApprox_2 *= knownTransmission(event, prevEvent, Ne - h);
						derivative += derivativeNeLogMulti(event, prevEvent);
						break;
					}
					derivativeApprox_1 *= biffurcation(event, prevEvent, Ne + h);
					derivativeApprox_2 *= biffurcation(event, prevEvent, Ne - h);
					derivative += derivativeNeLogBif(event, prevEvent);
					break;
				case MULTIFURCATION:
					// If event is recorded at time of transmission and
					// transmission tree branch is recipient,
					// log as transmission event.
					// This means no multiplication with rate \lambda*P_0
					if (!trNode.isRoot() && !donor && event.time == trNode.getParent().getHeight()) {
						derivativeApprox_1 *= knownTransmission(event, prevEvent, Ne + h);
						derivativeApprox_2 *= knownTransmission(event, prevEvent, Ne - h);
						derivative += derivativeNeLogMulti(event, prevEvent);
							break;
						}
					derivativeApprox_1 *= multifurcation(event, prevEvent, Ne + h);
					derivativeApprox_2 *= multifurcation(event, prevEvent, Ne - h);
					derivative += derivativeNeLogMulti(event, prevEvent);
				}
				prevEvent = event;
			}

			if (!trNode.isRoot() && prevEvent.time < trNode.getParent().getHeight()) {
						GeneTreeEvent mockEvent = new GeneTreeEvent();
				mockEvent.time = trNode.getParent().getHeight();
				mockEvent.lineages = prevEvent.lineages;
				derivativeApprox_1 *= interval(mockEvent, prevEvent, Ne + h);
				derivativeApprox_2 *= interval(mockEvent, prevEvent, Ne - h);
				derivative += derivativeNeLogInt(mockEvent, prevEvent);
				if (!donor) {
					derivativeApprox_1 *= knownTransmission(mockEvent, prevEvent, Ne + h);
					derivativeApprox_2 *= knownTransmission(mockEvent, prevEvent, Ne - h);
					derivative += derivativeNeLogMulti(mockEvent, prevEvent);
				}
			}
			double m = Math.log(derivativeApprox_1) - Math.log(derivativeApprox_2);
			m /= (2 * h);
			derivativesApprox.add(m);

			derivativesExact.add(derivative);

			if (derivativesExact.size() / (double) calcEvery == derivativesExact.size() / calcEvery
					&& derivativesExact.size() / calcEvery >= 1)
			{
				sumDerivativesApprox += derivativesApprox.subList((kl - 1) * calcEvery, kl * calcEvery).stream()
						.mapToDouble(a -> a)
						.sum();

				sumDerivativesExact += derivativesExact.subList((kl - 1) * calcEvery, kl * calcEvery).stream()
						.mapToDouble(a -> a)
						.sum();
				
				ps.print(derivativesExact.size() + "\t" + sumDerivativesApprox / (kl * calcEvery) + "\t"
						+ sumDerivativesExact / (kl * calcEvery));
				ps.print("\n");
				kl++;
			}
		}

		public void getAnalysis() {
			ps.close();
		}


		private double biffurcation(GeneTreeEvent event, GeneTreeEvent prevEvent, double Ne) {
			double ans = 0.0;
			ans += 1.0 / Ne;
			ans += (1.0 / waysToCoal(prevEvent.lineages, event.lineages))
					* gUp(prevEvent.lineages, event.lineages, tau, Ne)
					* lambda * P_0(event.time);


			return ans;
		}

		private double multifurcation(GeneTreeEvent event, GeneTreeEvent prevEvent, double Ne) {
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
		}

		private double knownTransmission(GeneTreeEvent event, GeneTreeEvent prevEvent, double Ne) {
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
		}

		private double interval(GeneTreeEvent event, GeneTreeEvent prevEvent, double Ne) {
			double ans = 0.0;
			ans -= prevEvent.lineages * (prevEvent.lineages - 1) * 0.5 * (1.0 / Ne) * (event.time - prevEvent.time);

			double sum = 0;
			for (int j = 1; j < prevEvent.lineages; j++) {
				sum += gUp(prevEvent.lineages, j, tau, Ne);

			}

			ans -= sum * lambda
					* integralP_0(prevEvent.time, event.time, lambda, mu, psi);

			return Math.exp(ans);
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

	/**
	 * Display error, print usage and exit with error.
	 */
	public static void printUsageAndError(String errMsg) {
		System.err.println(errMsg);
		System.err.println(helpMessage);
		System.exit(1);
	}

	/**
	 * Print usage info and exit.
	 */
	public static void printUsageAndExit() {
		System.out.println(helpMessage);
		System.exit(0);
	}


}
