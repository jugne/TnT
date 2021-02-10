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

	public Input<RealParameter> tauInput = new Input<>("tau",
			"Strength of the transmission bottleneck", Validate.REQUIRED);

	public Input<RealParameter> popSizesInput = new Input<RealParameter>("populationSizes",
			"Constant per-branch population sizes.", Validate.REQUIRED);
	public Input<RealParameter> originInput = new Input<RealParameter>("origin",
			"Start of branching process.", Validate.REQUIRED);

	// transformed parameters:
	public Input<RealParameter> expectedNInput = new Input<RealParameter>("expectedN",
			"The expected-N-at-present parameterisation of T", (RealParameter) null);
	public Input<RealParameter> diversificationRateInput = new Input<RealParameter>("diversificationRate",
			"Net diversification rate. Birth rate - death rate", Input.Validate.XOR, birthRateInput);
	public Input<Function> turnoverInput = new Input<Function>("turnover", "Turnover. Death rate/birth rate",
			Input.Validate.XOR, deathRateInput);
	public Input<RealParameter> samplingProportionInput = new Input<RealParameter>("samplingProportion",
			"The probability of sampling prior to death. Sampling rate/(sampling rate + death rate)",
			Input.Validate.XOR, samplingRateInput);

	// r parameter
	public Input<RealParameter> removalProbability = new Input<RealParameter>("removalProbability",
			"The probability that an individual is removed from the process after the sampling",
			(RealParameter) null);

	public Input<RealParameter> rhoProbability = new Input<RealParameter>("rho",
			"Probability of an individual to be sampled at present", (RealParameter) null);




	private double ploidy;

	private int transmissionNodeCount;

	private double lambda;
	private double tau;
	private double psi;
	private double mu;
	private double c1;
	private double c2;
	private double rho;
	private double r;
	protected boolean transform; // is true if the model is parametrised through transformed parameters
	RealParameter popSizes;
//	private double popSizes;
	private GeneTreeIntervals intervals;


	@Override
	public void initAndValidate() {
		if (birthRateInput.get() != null && deathRateInput.get() != null && samplingRateInput.get() != null) {

			transform = false;

		} else if (diversificationRateInput.get() != null && turnoverInput.get() != null
				&& samplingProportionInput.get() != null) {

			transform = true;

		} else {
			throw new IllegalArgumentException(
					"Either specify birthRate, deathRate and samplingRate OR specify diversificationRate, turnover and samplingProportion!");
		}

		ploidy = ploidyInput.get();
		intervals = geneTreeIntervalsInput.get();
		transmissionNodeCount = intervals.transmissionTree.getNodeCount();
//		popSizesInput.get().setDimension(transmissionNodeCount);
		popSizesInput.get().setDimension(2);
	}

	private void transformParameters() {
		double d = diversificationRateInput.get().getValue();
		double r_turnover = turnoverInput.get().getArrayValue();
		double s = samplingProportionInput.get().getValue();
		lambda = d / (1 - r_turnover);
		mu = r_turnover * lambda;
		psi = mu * s / (1 - s);
	}

	protected void updateParameters() {

		if (removalProbability.get() != null) {
			r = removalProbability.get().getValue();
		} else {
			r = 0.;
		}
		if (rhoProbability.get() != null) {
			rho = rhoProbability.get().getValue();
		} else {
			rho = 0.;
		}

		if (transform) {
			transformParameters();
		} else {
			lambda = birthRateInput.get().getValue();
			mu = deathRateInput.get().getArrayValue();
			psi = samplingRateInput.get().getValue();
		}

		c1 = Math.sqrt((lambda - mu - psi) * (lambda - mu - psi) + 4 * lambda * psi);
		c2 = -(lambda - mu - 2 * lambda * rho - psi) / c1;

		popSizes = popSizesInput.get();
		tau = tauInput.get().getValue();
	}

	@Override
	public double calculateLogP() {
//		lambda = birthRateInput.get().getValue();
//		mu = deathRateInput.get().getArrayValue();
//		psi = samplingRateInput.get().getValue();
//
//		if (rhoProbability.get() != null) {
//			rho = rhoProbability.get().getValue();
//		} else {
//			rho = 0.;
//		}
//		c1 = Math.sqrt((lambda - mu - psi) * (lambda - mu - psi) + 4 * lambda * psi);
//		c2 = -(lambda - mu - 2 * lambda * rho - psi) / c1;
//		
//		popSizes = popSizesInput.get();
//		tau = tauInput.get().getValue();

		updateParameters();


		logP = 0;

		// Calculate tree intervals
		HashMap<Integer, List<GeneTreeEvent>> eventList = intervals.getGeneTreeEventList();
		// if gene tree and species tree are incompatible
		if (eventList == null)
			return Double.NEGATIVE_INFINITY;
		SpeciesTreeInterface transmissionTree = (SpeciesTreeInterface) intervals.transmissionTreeInput.get()
				.getCurrent();

		
		for (Node trNode : transmissionTree.getNodesAsArray()) {
			boolean recipient = (!trNode.isRoot()
					&& trNode.getParent().getChild(0) != trNode && !trNode.getParent().isFake());
//			double popSize = popSizes;// popSizes.getValue(trNode.getNr());

			GeneTreeEvent prevEvent = new GeneTreeEvent();
			prevEvent.time = trNode.getHeight();

			for (GeneTreeEvent event : eventList.get(trNode.getNr())) {

				// Check that there are no events that leave non-positive number of lineages
				// this indicates that there is a bug in the package
				if (event.lineages == 0) {
					System.err.println(
							"Zero lineages after event!! Something is wrong with interval/event calculation. Please report this to the developers!");
					System.exit(1);
				}

				// Check if the event is at transmission time and on recipient side
				boolean eventAtTransmission = !trNode.isRoot() && recipient
						&& event.time == trNode.getParent().getHeight();

				// Contribution from every interval, except the last
				if (prevEvent.time < event.time) {
					logP += interval(event, prevEvent);
				}

				switch (event.type) {
				case SAMPLE:
					break;

				case BIFURCATION:
					// is bifurcation time is the end of transmission tree branch,
					// then it is the observed transmission event
					if (eventAtTransmission) {
						logP += transmission(event, prevEvent);
					} else {
						// otherwise, bifurcation contribution
						logP += biffurcation(event, prevEvent);
					}
					break;

				case MULTIFURCATION:
					// is multifurcation time is the end of transmission tree branch,
					// then it is the observed transmission event
					if (tau==0.0)
						logP = Double.NEGATIVE_INFINITY;
					if (eventAtTransmission) {
						logP += transmission(event, prevEvent);
					} else {
						// otherwise, multifurcation contribution
						logP += multifurcation(event, prevEvent);
					}
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
			// and transmission event at the end.
			if (logP != Double.NEGATIVE_INFINITY && !trNode.isRoot()
					&& prevEvent.time < trNode.getParent().getHeight()) {

				// mock event, to mark transmission that didn't result in merging of any
				// lineages
				GeneTreeEvent mockEvent = new GeneTreeEvent();
				mockEvent.time = trNode.getParent().getHeight();
				mockEvent.lineages = prevEvent.lineages;
				logP += interval(mockEvent, prevEvent);

				if (recipient) {
					logP += transmission(mockEvent, prevEvent);
				}
			}
		}

		return logP;
	}

	private double interval(GeneTreeEvent event, GeneTreeEvent prevEvent) {

		double ans = 0.0;


		// whole interval is after the origin, backwards in time
		if (prevEvent.time > originInput.get().getArrayValue(0)) {
			ans -= prevEvent.lineages * (prevEvent.lineages - 1) * 0.5 * (1.0 / (ploidy * popSizes.getValue(1)))
					* (event.time - prevEvent.time);
			return ans;
		}

		double end_time = event.time;
		if (event.time > originInput.get().getArrayValue(0)) {
			end_time = originInput.get().getArrayValue();
//			contribution after origin
			ans -= prevEvent.lineages * (prevEvent.lineages - 1) * 0.5 * (1.0 / (ploidy * popSizes.getValue(1)))
					* (event.time - end_time);
		}

//		contribution before origin
		ans -= prevEvent.lineages * (prevEvent.lineages - 1) * 0.5 * (1.0 / (ploidy * popSizes.getValue(0)))
				* (end_time - prevEvent.time);

		double sum = 0;
		sum = gUp(prevEvent.lineages, prevEvent.lineages, tau, popSizes.getValue(0));//

		ans -= (1 - sum) * lambda//
				* integralP_0(prevEvent.time, end_time);//

		return ans;
	}

	private double transmission(GeneTreeEvent event, GeneTreeEvent prevEvent) {
		double ans;
		double mult = 1.0;

		if (event.time > originInput.get().getArrayValue(0))
			return Double.NEGATIVE_INFINITY;

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
				* gUp(prevEvent.lineages, event.lineages, tau, popSizes.getValue(0));

		if (Math.log(ans) == Double.NEGATIVE_INFINITY)
			System.out.println();
		return Math.log(ans);
	}

	private double multifurcation(GeneTreeEvent event, GeneTreeEvent prevEvent) {
		double ans;
		double mult = 1.0;
		int sum = event.multiCoalSize.stream().mapToInt(Integer::intValue)
				.sum();

		int n_histories = event.multiCoalSize.size();

		if (event.time > originInput.get().getArrayValue(0) || tau == 0.0)
			return Double.NEGATIVE_INFINITY;
		for (int s = 0; s < n_histories; s++) {
			mult *= waysToCoal(event.multiCoalSize.get(s), 1);

			// W factor from NOAH A. ROSENBERG
			mult *= binomialInt(sum - (n_histories - s), event.multiCoalSize.get(s) - 1);
			sum -= event.multiCoalSize.get(s);
		}


		ans = (1.0 / waysToCoal(prevEvent.lineages, event.lineages)) * mult
				* gUp(prevEvent.lineages, event.lineages, tau, popSizes.getValue(0))
					* lambda * P_0(event.time);

		return Math.log(ans);
	}

	private double biffurcation(GeneTreeEvent event, GeneTreeEvent prevEvent) {
		double ans = 0.0;

		// event is after origin backwards in time
		if (event.time > originInput.get().getArrayValue(0)) {
			ans += 1.0 / (ploidy * popSizes.getValue(1));
			return Math.log(ans);
		}

		ans += 1.0 / (ploidy * popSizes.getValue(0));
		ans += (1.0 / waysToCoal(prevEvent.lineages, event.lineages))//
				* gUp(prevEvent.lineages, event.lineages, tau, popSizes.getValue(0))//
				* lambda * P_0(event.time);//

		return Math.log(ans);
	}


	public Node getRoot() {
		return intervals.geneTreeInput.get().getRoot();
	}


	private double gUp(int i, int j, double tau, double Ne) {
		double ans = 0.0;

		if (tau == 0.0 && i == j)
			return 1.0;
		if (tau == 0.0 && i != j)
			return 0.0;
		for (int k = j; k <= i; k++) {
			ans += (2.0 * k - 1) * Math.pow(-1.0, k - j) * f_1(j, k - 1) * f_2(i, k)
					* Math.exp(-(k * (k - 1) * tau * 0.5) / (ploidy * Ne)) /
					(MathUtils.factorialDouble(j) * MathUtils.factorialDouble(k - j) * f_1(i, k));
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

		return p0;
	}

	private double integralP_0(double t_0, double t_1) {

		double ans = (1.0 / (2.0 * lambda)) * ((t_1 - t_0) * (mu + psi + lambda - c1) + 2.0 * Math
				.log(((c2 - 1) * Math.exp(-c1 * t_0) - c2 - 1) / ((c2 - 1) * Math.exp(-c1 * t_1) - c2 - 1)));

		return ans;
	}

	/**
	 * @param i number of lineages before coalescent events, time increasing in the
	 *          past
	 * @param j number of lineages before after events, time increasing in the past
	 * @return nmber of ways i lineages can coalesce to j lineages
	 */
	private double waysToCoal(int i, int j) {
		if (i == j)
			return 1.0;
		double ans = MathUtils.factorialDouble(i);
		ans *= MathUtils.factorialDouble(i - 1);
		ans /= Math.pow(2, i - j);
		ans /= MathUtils.factorialDouble(j);
		ans /= MathUtils.factorialDouble(j - 1);

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
