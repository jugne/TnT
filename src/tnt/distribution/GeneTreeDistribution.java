package tnt.distribution;

import java.util.HashMap;
import java.util.List;
import java.util.Random;

import org.apache.commons.math.util.MathUtils;

import beast.core.Distribution;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import starbeast2.SpeciesTreeInterface;

public class GeneTreeDistribution extends Distribution {

	public Input<GeneTreeIntervals> geneTreeIntervalsInput = new Input<>("geneTreeIntervals",
			"intervals for a gene tree", Validate.REQUIRED);
	public Input<Double> ploidyInput = new Input<>("ploidy",
			"Ploidy (copy number) for this gene, typically a whole number or half (default is 1).", 1.0);
	public Input<RealParameter> birthRateInput = new Input<RealParameter>("birthRate", "Birth rate",
			Input.Validate.REQUIRED);
	public Input<Function> deathRateInput = new Input<Function>("deathRate", "Death rate", Input.Validate.REQUIRED);
	public Input<RealParameter> samplingRateInput = new Input<RealParameter>("samplingRate",
			"Sampling rate per individual", Input.Validate.REQUIRED);

	// r parameter
	public Input<RealParameter> removalProbability = new Input<RealParameter>("removalProbability",
			"The probability that an individual is removed from the process after the sampling",
			Input.Validate.REQUIRED);

	public Input<RealParameter> rhoProbability = new Input<RealParameter>("rho",
			"Probability of an individual to be sampled at present", (RealParameter) null);

	public Input<RealParameter> tauInput = new Input<>("tau",
			"Strength of the transmission bottleneck", Validate.REQUIRED);

	public Input<RealParameter> popSizesInput = new Input<RealParameter>("populationSizes",
			"Constant per-branch population sizes.", Validate.REQUIRED);


	private double ploidy;

	private int transmissionNodeCount;

	private double lambda;
	private double tau;
	private double psi;
	private double mu;
	private double c1;
	private double c2;
	private double r;
	private double rho;
	RealParameter popSizes;
	private GeneTreeIntervals intervals;


	@Override
	public void initAndValidate() {
		ploidy = ploidyInput.get();
		intervals = geneTreeIntervalsInput.get();
		transmissionNodeCount = intervals.transmissionTree.getNodeCount();
		popSizesInput.get().setDimension(transmissionNodeCount);
	}

	@Override
	public double calculateLogP() {
		lambda = birthRateInput.get().getValue();
		mu = deathRateInput.get().getArrayValue();
		psi = samplingRateInput.get().getValue();

		r = removalProbability.get().getValue();
		if (rhoProbability.get() != null) {
			rho = rhoProbability.get().getValue();
		} else {
			rho = 0.;
		}
		c1 = Math.sqrt((lambda - mu - psi) * (lambda - mu - psi) + 4 * lambda * psi);
		c2 = -(lambda - mu - 2 * lambda * rho - psi) / c1;

		popSizes = popSizesInput.get();
		tau = tauInput.get().getValue();

		logP = 0;

		// Calculate tree intervals
		HashMap<Integer, List<GeneTreeEvent>> eventList = intervals.getGeneTreeEventList();
		SpeciesTreeInterface transmissionTree = (SpeciesTreeInterface) intervals.transmissionTreeInput.get()
				.getCurrent();

		// TODO check if this is working as intended
		// if gene tree and species tree are incompatible
		if (eventList == null)
			return Double.NEGATIVE_INFINITY;



		
		for (Node trNode : transmissionTree.getNodesAsArray()) {
//			if (trNode.isDirectAncestor())
//				continue;
			// define donor as a left child
//			boolean donor = (trNode.isRoot()
//					|| (trNode.getParent().getChild(0) == trNode && !trNode.getParent().isFake()));
			boolean recipient = (!trNode.isRoot()
					&& trNode.getParent().getChild(0) != trNode && !trNode.getParent().isFake());
			double popSize = popSizes.getValue(trNode.getNr());

			GeneTreeEvent prevEvent = new GeneTreeEvent();
			prevEvent.time = trNode.getHeight();

			for (GeneTreeEvent event : eventList.get(trNode.getNr())) {

				// Check that there are no events that leave non-positive number of lineages
				if (event.lineages == 0) {
					System.err
							.println("Zero lineages after event!! Something is wrong with interval/event calculation");
					System.exit(1);
				}

				// Contribution from every interval, except the last
				if (prevEvent.time < event.time) {
					logP += interval(event, prevEvent, popSize);
				}

				switch (event.type) {
				case SAMPLE:
					break;

				case BIFURCATION:
					// is bifurcation time is the end of transmission tree branch,
					// then it is the observed transmission event
					if (!trNode.isRoot() && recipient && event.time == trNode.getParent().getHeight()) {
						logP += transmission(event, prevEvent, popSize);
					}
					// otherwise, bifurcation contribution
					logP += biffurcation(event, prevEvent, popSize);
					break;

				case MULTIFURCATION:
					// is multifurcation time is the end of transmission tree branch,
					// then it is the observed transmission event
					if (!trNode.isRoot() && recipient && event.time == trNode.getParent().getHeight()) {
						logP += transmission(event, prevEvent, popSize);
					}
					logP += multifurcation(event, prevEvent, popSize);
					break;


				}

				if (logP == Double.NEGATIVE_INFINITY)
					break;

				// set previous event to the current
				prevEvent = event;
			}

			// if after going through all the events, there is still an interval left
			// between last event
			// on a gene tree and transmission tree branch, we have to calculate
			// contribution of this interval
			// and transmission venet at the end.
			if (logP != Double.NEGATIVE_INFINITY && !trNode.isRoot()
					&& prevEvent.time < trNode.getParent().getHeight()) {

				// mock event, to mark transmission that didn't result in merging of any
				// lineages
				GeneTreeEvent mockEvent = new GeneTreeEvent();
				mockEvent.time = trNode.getParent().getHeight();
				mockEvent.lineages = prevEvent.lineages;
				logP += interval(mockEvent, prevEvent, popSize);

				if (recipient) {
					logP += transmission(mockEvent, prevEvent, popSize);
				}
			}
		}

		return logP;
	}

	private double interval(GeneTreeEvent event, GeneTreeEvent prevEvent, double Ne) {

		double ans = 0.0;
		ans -= prevEvent.lineages * (prevEvent.lineages - 1) * 0.5 * (1.0 / (ploidy * Ne))
				* (event.time - prevEvent.time);

		double sum = 0;
		sum = gUp(prevEvent.lineages, prevEvent.lineages, tau, Ne);

		ans -= (1 - sum) * lambda
				* integralP_0(prevEvent.time, event.time);

		return ans;
	}

	private double transmission(GeneTreeEvent event, GeneTreeEvent prevEvent, double Ne) {
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

		return Math.log(ans);
	}

	private double multifurcation(GeneTreeEvent event, GeneTreeEvent prevEvent, double Ne) {
		double ans;
		double mult = 1.0;
		int sum = event.multiCoalSize.stream().mapToInt(Integer::intValue)
				.sum();
		int n_histories = event.multiCoalSize.size();

		for (int s = 0; s < n_histories; s++) {
			mult *= waysToCoal(event.multiCoalSize.get(s), 1);
			// W factor from NOAH A. ROSENBERG
			mult *= binomialInt(sum - (n_histories - s), event.multiCoalSize.get(s) - 1);
			sum -= event.multiCoalSize.get(s);
		}
		ans = (1.0 / waysToCoal(prevEvent.lineages, event.lineages)) * mult
				* gUp(prevEvent.lineages, event.lineages, tau, Ne)
				* lambda * P_0(event.time);

		return Math.log(ans);
	}

	private double biffurcation(GeneTreeEvent event, GeneTreeEvent prevEvent, double Ne) {
		double ans = 0.0;
		ans += 1.0 / (ploidy * Ne);
		ans += (1.0 / waysToCoal(prevEvent.lineages, event.lineages))
				* gUp(prevEvent.lineages, event.lineages, tau, Ne)
				* lambda * P_0(event.time);


		return Math.log(ans);
	}


	public Node getRoot() {
		return intervals.geneTreeInput.get().getRoot();
	}


	private double gUp(int i, int j, double tau, double Ne) {
		double ans = 0.0;
//		if (tau == 0)
//			return ans;
		for (int k = j; k <= i; k++) {
			ans += (2 * k - 1) * Math.pow(-1, k - j) * f_1(j, k - 1) * f_2(i, k)
					* Math.exp(-(k * (k - 1) * tau * 0.5) / (ploidy * Ne)) /
					(MathUtils.factorial(j) * MathUtils.factorial(k - j) * f_1(i, k));
		}

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

	private double P_0(double t) {

		double p0 = (lambda + mu + psi
				+ c1 * ((Math.exp(-c1 * t) * (1 - c2) - (1 + c2))
						/ (Math.exp(-c1 * t) * (1 - c2) + (1 + c2))))
				/ (2.0 * lambda);

		return r + (1 - r) * p0;
	}

	private double integralP_0(double t_0, double t_1) {

		double ans = (1.0 / (2.0 * lambda)) * ((t_1 - t_0) * (mu + psi + lambda - c1) + 2.0 * Math
				.log(((c2 - 1) * Math.exp(-c1 * t_0) - c2 - 1) / ((c2 - 1) * Math.exp(-c1 * t_1) - c2 - 1)));

		return r + (1 - r) * ans;
	}

	private double waysToCoal(int i, int j) {
		double ans = MathUtils.factorial(i);
		ans *= MathUtils.factorial(i - 1);
		ans /= Math.pow(2, i - j);
		ans /= MathUtils.factorial(j);
		ans /= MathUtils.factorial(j - 1);

		return ans;
	}

	private long binomialInt(int n, int k) {
		if (k == n)
			return 1;
		if (k > n - k)
			k = n - k;

		long binom = 1;
		for (int i = 1; i <= k; i++)
			binom = binom * (n + 1 - i) / i;
		return binom;
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
