package tnt.distribution;

import java.util.HashMap;
import java.util.List;
import java.util.Random;

import org.apache.commons.math.util.MathUtils;

import bdmmprime.parameterization.Parameterization;
import beast.core.Distribution;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import starbeast2.SpeciesTreeInterface;
import tnt.util.Tools;

public class SkyGeneTreeDistribution extends Distribution {

	public Input<GeneTreeIntervals> geneTreeIntervalsInput = new Input<>("geneTreeIntervals",
			"intervals for a gene tree", Validate.REQUIRED);
	public Input<Double> ploidyInput = new Input<>("ploidy",
			"Ploidy (copy number) for this gene, typically a whole number or half (default is 1).", 1.0);

	public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
			"BDMM parameterization",
			Input.Validate.REQUIRED);

	public Input<RealParameter> pairwiseProbBottleneckInput = new Input<>("pairwiseProbBottleneck",
			"Pairwise coalescent probability in the bottleneck");

	public Input<RealParameter> durationInput = new Input<>("bottleneckDuration",
			"Duration of the transmission bottleneck. Use either duration or strength (tau).");

	public Input<RealParameter> tauInput = new Input<>("tau",
			"Strength of the transmission bottleneck");

	public Input<RealParameter> popSizesInput = new Input<RealParameter>("populationSizes",
			"Constant per-branch population sizes.", Validate.REQUIRED);


	// transformed parameters:
	public Input<RealParameter> expectedNInput = new Input<RealParameter>("expectedN",
			"The expected-N-at-present parameterisation of T", (RealParameter) null);
	
    public Input<Function> finalSampleOffsetInput = new Input<>("finalSampleOffset",
            "If provided, the difference in time between the final sample and the end of the BD process.",
            new RealParameter("0.0"));

	public Input<Boolean> popSizePerBranchInput = new Input<>("popSizePerBranch",
			"If pop sizes are estimated for each transmission tree branch. "
					+ "Default false: only before and after origin sizes are estimated.",
			false);

	private double ploidy;

	private double bottleneckDuration;

	private Parameterization parameterization;

	private double rho_i;
	private double lambda_i;
	private double mu_i;
	private double psi_i;
	private double t_i;

	private double[] A;
	private double[] B;
	private double A_i;
	private double B_i;
	private double origin;
	private Function finalSampleOffset;

	private double popSizePerBranch;
	private double popSizeOrigin;

	private int popSizeDim;

	SpeciesTreeInterface transmissionTree;

	private GeneTreeIntervals intervals;


	@Override
	public void initAndValidate() {

		ploidy = ploidyInput.get();
		intervals = geneTreeIntervalsInput.get();

		transmissionTree = intervals.transmissionTreeInput.get();

		popSizeDim = popSizePerBranchInput.get() ? transmissionTree.getNodeCount() + 1 : 2;
		popSizesInput.get().setDimension(popSizeDim);

		parameterization = parameterizationInput.get();
		
		finalSampleOffset = finalSampleOffsetInput.get();

		A = new double[parameterization.getTotalIntervalCount()];
		B = new double[parameterization.getTotalIntervalCount()];
	}


	protected void updateParameters(int trNodeNr) {
		
		popSizePerBranch = popSizePerBranchInput.get() ? popSizesInput.get().getValue(trNodeNr)
				: popSizesInput.get().getValue(0);

		if (tauInput.get() != null)
			bottleneckDuration = tauInput.get().getValue() * ploidy * popSizePerBranch;
		else if (durationInput.get() != null)
			bottleneckDuration = durationInput.get().getValue();
		else if (pairwiseProbBottleneckInput.get() != null) {
			if (pairwiseProbBottleneckInput.get().getValue() == 0)
				bottleneckDuration = 0;
			else if (pairwiseProbBottleneckInput.get().getValue() > 1.0
					|| pairwiseProbBottleneckInput.get().getValue() < 0)
				throw new Error("pairwiseProbBottleneck must be between 0 and 1.0");
			else {
				bottleneckDuration = -Math.log(1 - pairwiseProbBottleneckInput.get().getValue()) * ploidy
						* popSizePerBranch;
			}
		} else
			throw new Error("Either bottleneckDuration, tau or pairwiseProbBottleneck must be specified!");

		popSizeOrigin = popSizesInput.get().getValue(popSizeDim-1);
	}

	private void computeConstants(double[] A, double[] B) {

		for (int i = parameterization.getTotalIntervalCount() - 1; i >= 0; i--) {

			double p_i_prev;
			if (i + 1 < parameterization.getTotalIntervalCount()) {
				p_i_prev = get_p_i(parameterization.getBirthRates()[i + 1][0],
						parameterization.getDeathRates()[i + 1][0],
						parameterization.getSamplingRates()[i + 1][0],
						A[i + 1], B[i + 1],
						parameterization.getIntervalEndTimes()[i + 1],
						parameterization.getIntervalEndTimes()[i]);
			} else {
				p_i_prev = 1.0;
			}

			double rho_i = parameterization.getRhoValues()[i][0];
			double lambda_i = parameterization.getBirthRates()[i][0];
			double mu_i = parameterization.getDeathRates()[i][0];
			double psi_i = parameterization.getSamplingRates()[i][0];

			A[i] = Math.sqrt((lambda_i - mu_i - psi_i) * (lambda_i - mu_i - psi_i) + 4 * lambda_i * psi_i);
			B[i] = ((1 - 2 * (1 - rho_i) * p_i_prev) * lambda_i + mu_i + psi_i) / A[i];
		}
	}


	@Override
	public double calculateLogP() {

		logP = 0;

		// Calculate tree intervals
		HashMap<Integer, List<GeneTreeEvent>> eventList = intervals.getGeneTreeEventList();
		// if gene tree and species tree are incompatible
		if (eventList == null)
			return Double.NEGATIVE_INFINITY;


		origin = parameterization.originInput.get().getArrayValue(0);

		computeConstants(A, B);

		
		for (Node trNode : transmissionTree.getNodesAsArray()) {
			if (eventList.get(trNode.getNr()) != null) {
			updateParameters(trNode.getNr());

			int i = parameterization.getIntervalIndex(parameterization.getNodeTime(trNode, finalSampleOffset.getArrayValue()));
			double int_end_time = Double.POSITIVE_INFINITY;
			if (i != 0)
				int_end_time = origin - parameterization.getIntervalEndTimes()[i - 1];

			rho_i = parameterization.getRhoValues()[i][0];
			lambda_i = parameterization.getBirthRates()[i][0];
			mu_i = parameterization.getDeathRates()[i][0];
			psi_i = parameterization.getSamplingRates()[i][0];
				t_i = parameterization.getIntervalEndTimes()[i];
				A_i = A[i];
				B_i = B[i];


			boolean recipient = (!trNode.isRoot()
					&& trNode.getParent().getChild(0) != trNode && !trNode.getParent().isFake());

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
							&& Tools.equalWithPrecisionDouble(event.time, trNode.getParent().getHeight());

				// Contribution from every interval, except the last				
				if (prevEvent.time < event.time) {
					double startTime = prevEvent.time;
					if (event.time > int_end_time) {
						while (event.time > int_end_time && event.time < origin) {
							logP += interval(int_end_time, startTime, prevEvent);
							startTime = int_end_time;
							i = parameterization.getIntervalIndex(origin - startTime - 0.00001);
							int_end_time = Double.POSITIVE_INFINITY;
							if (i != 0)
								int_end_time = origin - parameterization.getIntervalEndTimes()[i - 1];
							rho_i = parameterization.getRhoValues()[i][0];
							lambda_i = parameterization.getBirthRates()[i][0];
							mu_i = parameterization.getDeathRates()[i][0];
							psi_i = parameterization.getSamplingRates()[i][0];
								t_i = parameterization.getIntervalEndTimes()[i];
								A_i = A[i];
								B_i = B[i];
						}
						if (event.time != startTime) {
							logP += interval(event.time, startTime, prevEvent);
						}
					} else {
						logP += interval(event.time, startTime, prevEvent);
					}
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
					if (bottleneckDuration == 0.0)
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
				double startTime = prevEvent.time;
				if (mockEvent.time > int_end_time) {
					while (mockEvent.time > int_end_time) {
						logP += interval(int_end_time, startTime, prevEvent);
						startTime = int_end_time;
						i = parameterization.getIntervalIndex(origin - startTime - 0.00001);
						int_end_time = Double.POSITIVE_INFINITY;
						if (i != 0)
							int_end_time = origin - parameterization.getIntervalEndTimes()[i - 1];
						rho_i = parameterization.getRhoValues()[i][0];
						lambda_i = parameterization.getBirthRates()[i][0];
						mu_i = parameterization.getDeathRates()[i][0];
						psi_i = parameterization.getSamplingRates()[i][0];
							t_i = parameterization.getIntervalEndTimes()[i];
							A_i = A[i];
							B_i = B[i];
					}
					if (mockEvent.time != startTime) {
						logP += interval(mockEvent.time, startTime, prevEvent);
					}
				} else {
					logP += interval(mockEvent.time, startTime, prevEvent);
				}

				if (recipient) {
					logP += transmission(mockEvent, prevEvent);
				}
			}
		}
		}

		return logP;
	}

	private double interval(double endTime, double startTime, GeneTreeEvent prevEvent) {

		double ans = 0.0;


		// whole interval is after the origin, backwards in time
		if (startTime > origin) {
			ans -= prevEvent.lineages * (prevEvent.lineages - 1) * 0.5 * (1.0 / (ploidy * popSizeOrigin))
					* (endTime - startTime);
			return ans;
		}

		double end_time = endTime;
		if (endTime > origin) {
			end_time = origin;
//			contribution after origin
			ans -= prevEvent.lineages * (prevEvent.lineages - 1) * 0.5 * (1.0 / (ploidy * popSizeOrigin))
					* (endTime - end_time);
		}

//		contribution before origin
		ans -= prevEvent.lineages * (prevEvent.lineages - 1) * 0.5 * (1.0 / (ploidy * popSizePerBranch))
				* (end_time - prevEvent.time);

		double sum = 0;
		sum = gUp(prevEvent.lineages, prevEvent.lineages, bottleneckDuration, popSizePerBranch);//

		ans -= (1 - sum) * lambda_i//
				* integral_p_i(parameterization.getAge(startTime, finalSampleOffset.getArrayValue()),
						parameterization.getAge(end_time, finalSampleOffset.getArrayValue()));//

		return ans;
	}

	private double transmission(GeneTreeEvent event, GeneTreeEvent prevEvent) {
		double ans;
		double mult = 1.0;

		if (event.time > origin)
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
				* gUp(prevEvent.lineages, event.lineages, bottleneckDuration, popSizePerBranch);

//		if (Math.log(ans) == Double.NEGATIVE_INFINITY)
//			System.out.println();
		return Math.log(ans);
	}

	private double multifurcation(GeneTreeEvent event, GeneTreeEvent prevEvent) {
		double ans;
		double mult = 1.0;
		int sum = event.multiCoalSize.stream().mapToInt(Integer::intValue)
				.sum();

		int n_histories = event.multiCoalSize.size();

		if (event.time > origin || bottleneckDuration == 0.0)
			return Double.NEGATIVE_INFINITY;
		for (int s = 0; s < n_histories; s++) {
			mult *= waysToCoal(event.multiCoalSize.get(s), 1);

			// W factor from NOAH A. ROSENBERG
			mult *= binomialInt(sum - (n_histories - s), event.multiCoalSize.get(s) - 1);
			sum -= event.multiCoalSize.get(s);
		}


		ans = (1.0 / waysToCoal(prevEvent.lineages, event.lineages)) * mult
				* gUp(prevEvent.lineages, event.lineages, bottleneckDuration, popSizePerBranch)
				* lambda_i * get_p_i(lambda_i, mu_i, psi_i, A_i, B_i, t_i,
						parameterization.getAge(event.time, finalSampleOffset.getArrayValue()));

		return Math.log(ans);
	}

	private double biffurcation(GeneTreeEvent event, GeneTreeEvent prevEvent) {
		double ans = 0.0;

		// event is after origin backwards in time
		if (event.time > origin) {
			ans += 1.0 / (ploidy * popSizeOrigin);
			return Math.log(ans);
		}

		ans += 1.0 / (ploidy * popSizePerBranch);
		ans += (1.0 / waysToCoal(prevEvent.lineages, event.lineages))//
				* gUp(prevEvent.lineages, event.lineages, bottleneckDuration, popSizePerBranch)//
				* lambda_i * get_p_i(lambda_i, mu_i, psi_i, A_i, B_i, t_i,
						parameterization.getAge(event.time, finalSampleOffset.getArrayValue()));//

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
			ans += (2.0 * k - 1) * Math.pow(-1.0, k - j) * Math.exp(-(k * (k - 1) * tau * 0.5) / (ploidy * Ne))
					* (f_1(j, k - 1) / MathUtils.factorialDouble(k - j)) * (f_2(i, k) / MathUtils.factorialDouble(j))
					/
					f_1(i, k);
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
	
	private double get_p_i(double lambda, double mu, double psi, double A, double B, double t_i, double t) {

		if (lambda > 0.0) {
			double v = Math.exp(A * (t_i - t)) * (1 + B);
			double ans = (lambda + mu + psi - A * (v - (1 - B)) / (v + (1 - B)))
					/ (2 * lambda);
			return ans;
        } else {
            // The limit of p_i as lambda -> 0
            return 0.5;
        }
    }

	private double integral_p_i(double t_0, double t_1) {
		double t0 = t_i - t_0;
		double t1 = t_i - t_1;
		double ans = (1.0 / (2.0 * lambda_i)) * ((t1 - t0) * (mu_i + psi_i + lambda_i + A_i) + 2.0 * Math
				.log(((-B_i - 1) * Math.exp(A_i * t0) + B_i - 1)
						/ ((-B_i - 1) * Math.exp(A_i * t1) + B_i - 1)));

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
