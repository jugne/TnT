package tnt.likelihood;

import java.util.ArrayList;
import java.util.List;

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
//		Taxon t17_1 = new Taxon();
//		t17_1.setID("t17_1");
//		Taxon t18_1 = new Taxon();
//		t18_1.setID("t18_1");
//		Taxon t20_1 = new Taxon();
//		t20_1.setID("t20_1");
//		Taxon t29_1 = new Taxon();
//		t29_1.setID("t29_1");
//		Taxon t33_1 = new Taxon();
//		t33_1.setID("t33_1");

		taxonSet.initByName("taxon", t7_1, "taxon", t10_1, "taxon", t11_1, "taxon", t2_1, "taxon", t22_1, "taxon",
				t6_1);// , "taxon", t17_1, "taxon", t18_1, "taxon", t20_1, "taxon", t29_1, "taxon",
						// t33_1);

		TreeParser transmissionTreeInput = new TreeParser();
		transmissionTreeInput.initByName("newick", speciesTreeNewick, "adjustTipHeights", false, "IsLabelledNewick",
				true);

		TraitSet sampleCounts = new TraitSet();
		sampleCounts.initByName("traitname", "sampleCounts", "taxa", taxonSet, "value", " t7_1=2,\n" +
				"t10_1=5,\n" +
				"t11_1=4,\n" +
				"t6_1=2,\n" +
				"t2_1=6,\n" +
				"t22_1=5");

		RealParameter populationSizes = new RealParameter("1.0");
		RealParameter bottleneckStrength = new RealParameter("0.0");

		RealParameter birthRate = new RealParameter("1.0");
		RealParameter deathRate = new RealParameter("0.4");
		RealParameter samplingRate = new RealParameter("0.3");

		SimulatedGeneTree geneTree = new SimulatedGeneTree();
		geneTree.setID("gene_tree_truth");
		geneTree.initByName("complete", true, "transmissionTreeInput", transmissionTreeInput,
				"sampleCounts", sampleCounts,
				"populationSizes", populationSizes,
				"bottleneckStrength", bottleneckStrength,
				"birthRate", birthRate,
				"deathRate", deathRate,
				"samplingRate", samplingRate);

		GeneTreeIntervals intervals = new GeneTreeIntervals();
		intervals.initByName("geneTree", geneTree);

		Integer runs = 1000;
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

		List<Double> derivatives = new ArrayList<Double>();

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

//			intervals.initAndValidate();

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
			List<GeneTreeEvent> eventList = intervals.getGeneTreeEventList();
			
			double derivative = 0;
//			double dt = 0.00001;
//			double time = 0.0;
			GeneTreeEvent prevEvent = new GeneTreeEvent();
			prevEvent.time = 0.0;
			
			for (int i=0; i < eventList.size(); i++) {
				GeneTreeEvent event = eventList.get(i);
				if (prevEvent.time < event.time) {
					derivative += derivativeNeInterval(event, prevEvent);
				}
				switch (event.type) {
				case SAMPLE:
					break;
				case BIFURCATION:
					derivative += derivativeNeBifurcation(event);
					break;
				case MULTIFURCATION:
					derivative += derivativeNeMultifurcation(event, prevEvent);
				}
				prevEvent = event;

			}

			derivatives.add(derivative);
		}

		public List<Double> getAnalysis() {
			return derivatives;
		}

		private double derivativeNeMultifurcation(GeneTreeEvent event, GeneTreeEvent prevEvent) {
			
			int i = prevEvent.lineages;
			int j = event.lineages;
			
			double ans = 0.0;
			ans += derivativeNeG(i, j, tau, Ne);
			ans *= lambda;
			ans *= P_0(event.time, lambda, mu, psi);
			ans *= (1.0 / waysToCoal(i, j));
					
			
			return 0.0;
		}

		private double derivativeNeG(int i, int j, double tau, double Ne) {
			
			double sum = 0;
			double sumDerivatives = 0;
			for (int s = 1; s <= i; s++) {
				sum += gUp(i, s, tau, Ne);
				sumDerivatives += derivativeNeGUp(i, s, tau, Ne);
			}

			double ans = (derivativeNeGUp(i, j, tau, Ne) * sum);
			ans -= (gUp(i, j, tau, Ne) * sumDerivatives);
			ans /= Math.pow(sum, 2);
			
			return ans;
		}

		private double derivativeNeGUp(int i, int j, double tau, double Ne) {

			double ans = 0.0;
			for (int k = j; k <= i; k++) {
				ans += (2 * k - 1) * Math.pow(-1, k - j) * f_1(j, k - 1) * f_2(i, k)
						* Math.exp(-(k * (k - 1) * tau * 0.5) / Ne) * ((k * (k - 1) * tau * 0.5) / Math.pow(Ne, 2)) /
						(MathUtils.factorial(j) * MathUtils.factorial(k - j) * f_1(i, k));
			}

			return ans;
		}

		private double derivativeNeInterval(GeneTreeEvent event, GeneTreeEvent prevEvent) {
			
//			double sum = 0;
//			for (int j = 1; j < prevEvent.lineages; j++) {
//				sum += derivativeNeG(prevEvent.lineages, j, tau, Ne) * waysToCoal(prevEvent.lineages, j);
//			}
			
			return (prevEvent.lineages) * (prevEvent.lineages - 1) * 0.5 * (event.time - prevEvent.time)
					* (1.0 / Math.pow(Ne, 2));// - (sum * lambda * integralP_0(prevEvent.time, event.time, lambda, mu,
											// psi));
		}

		private double derivativeNeBifurcation(GeneTreeEvent event) {

			return -1.0 / Ne; // + derivativeNeG(event.lineages + 1, event.lineages, tau, Ne)* lambda
//					* P_0(event.time, lambda, mu, psi)) *
//					(1.0 / ((1.0 / Ne) + g(event.lineages + 1, event.lineages, tau, Ne) * lambda
//							* P_0(event.time, lambda, mu, psi)));

		}
		
		private double g(int i, int j, double tau, double Ne) {
			double normalised = 0.0;
			
			for (int s = 1; s <= i; s++) {
				normalised += gUp(i, s, tau, Ne);
			}

			double ans = gUp(i, j, tau, Ne);
			ans /= normalised;

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

		private double P_0(double t, double lambda, double mu, double psi) {

			double c_1 = Math.abs(Math.sqrt(Math.pow((lambda - mu - psi), 2) + 4 * lambda * psi));
			double c_2 = -(lambda - mu - psi) / c_1;

			
			return (lambda + mu + psi
					+ c_1 * ((Math.exp(c_1 * t) * (1 - c_2) - (1 + c_2)) / (Math.exp(c_1 * t) * (1 - c_2) + (1 + c_2))))
					/ (2.0 * lambda);
		}

		private double integralP_0(double t_0, double t_1, double lambda, double mu, double psi) {
			double c_1 = Math.abs(Math.sqrt(Math.pow((lambda - mu - psi), 2) + 4 * lambda * psi));
			double c_2 = -(lambda - mu - psi) / c_1;

			double ans = (1.0 / (2.0 * lambda)) * ((t_1 - t_0) * (mu + psi + lambda - c_1) + 2.0 * Math
					.log(((c_2 - 1) * Math.exp(-c_1 * t_0) - c_2 - 1) / (c_2 - 1) * Math.exp(c_1 * t_1) - c_2 - 1));

			return ans;
		}

		private double waysToCoal(int i, int j) {
			double ans = MathUtils.factorial(i);
			ans *= MathUtils.factorial(i - 1);
			ans /= Math.pow(2, i - j);
			ans /= MathUtils.factorial(j);
			ans /= MathUtils.factorial(j - 1);

			return 1.0;
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

	}

}
