package tnt.mapping;


import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import bdmmprime.parameterization.Parameterization;
import beast.core.BEASTObject;
import beast.core.Function;
import beast.core.Input;
import beast.core.StateNode;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.Parameter;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import starbeast2.PopulationModel;
import tnt.logger.TransmissionTreeLogger;
import tnt.transmissionTree.TransmissionTree;


/**
 * Maps hidden host switch events on a birth-death transmission tree.
 * 
 * @author Ugne Stolz <ugne.stolz@protonmail.com>
 * @date 7 Mar 2022
 */
public class MappedTree extends TransmissionTree {

	public Input<Boolean> mapOnInitInput = new Input<>("mapOnInit",
			"If true, mapping will be performed when object is " +
					"first initialize.",
			true);

	public Input<TransmissionTree> treeInput = new Input<>("tree", "Tree with no hidden host switch events.",
			Input.Validate.REQUIRED);

	public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
			"BDMM parameterization",
			Input.Validate.REQUIRED);

	// transformed parameters:
	public Input<Function> finalSampleOffsetInput = new Input<>("finalSampleOffset",
			"If provided, the difference in time between the final sample and the end of the BD process.",
			new RealParameter("0.0"));

	public Input<Boolean> remapOnLogInput = new Input<>("remapOnLog",
			"If true, mapping will be regenerated when this object " +
					"is logged.",
			true);
	
	public Input<Function> hiddenEventsCounterInput = new Input<>("nHiddenEvents",
			"Number of hidden events. Set by sampler.",
			new IntegerParameter("0"));

	public Input<List<TaxonSet>> taxonsetsInput = new Input<>("taxonsets",
			"a separate list of taxa for samples collected from the same patient", new ArrayList<>());





	private TransmissionTree unmappedTree;

	private Parameterization parameterization;

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

	private int nHiddenEvents;

	private List<String> constrainedTaxons = new ArrayList<>();
	private DecimalFormat df;



	@Override
	public void initAndValidate() {
		treeInput.get().initAndValidate();
		parameterization = parameterizationInput.get();
		origin = parameterization.originInput.get().getArrayValue(0);
		finalSampleOffset = finalSampleOffsetInput.get();
		unmappedTree = treeInput.get();

		for (TaxonSet t : taxonsetsInput.get()) {
			constrainedTaxons.addAll(t.asStringList());
		}

		A = new double[parameterization.getTotalIntervalCount()];
		B = new double[parameterization.getTotalIntervalCount()];

		if (mapOnInitInput.get())
			map();
	}

	public void map() {
		computeConstants(A, B);
		origin = parameterization.originInput.get().getArrayValue(0);
		unmappedTree = treeInput.get();
		unmappedTree.orientateTree();
		nHiddenEvents = 0;

		boolean skip = false;
		Node typedRoot = paintHiddenEvents(unmappedTree.getRoot(), skip);

		// Ensure internal nodes are numbered correctly. (Leaf node numbers and
		// labels are matched to those in the untyped tree during the simulation.)
		numberInternalNodesOnSubtree(typedRoot, unmappedTree.getLeafNodeCount());

		assignFromWithoutID(new Tree(typedRoot));

		IntegerParameter hiddenEventsCounter = (IntegerParameter) hiddenEventsCounterInput.get();
		hiddenEventsCounter.setValue(nHiddenEvents);

	}

	private Node paintHiddenEvents(Node subroot, boolean skip) {

		Node currentNode = new Node();

		if(!skip && subroot.isFake() && constrainedTaxons.contains(subroot.getDirectAncestorChild().getID())){
			skip=true;
		}
		double t_end_branch = origin - finalSampleOffset.getArrayValue();

		if (!subroot.isRoot())
			t_end_branch = subroot.getParent().getHeight();
		double t_start_branch = subroot.getHeight();

		int i = parameterization.getNodeIntervalIndex(subroot, finalSampleOffset.getArrayValue());
		double t_end_int = origin - finalSampleOffset.getArrayValue();
		if (i != 0)
			t_end_int = origin - parameterization.getIntervalEndTimes()[i - 1];

		updateParameters(i);
		
		double t_start = t_start_branch;
		List<Double> sampledTimes = new ArrayList<>();
		if (!skip) {
			while (t_end_branch > t_end_int) {
				// record event times
				sampledTimes.addAll(sampleTimes(t_start, t_end_int));

				i = i - 1;
				t_start = t_end_int;
				t_end_int = parameterization.getIntervalEndTimes()[i - 1];
				updateParameters(i);
			}

			// record event times
			sampledTimes.addAll(sampleTimes(t_start, t_end_branch));
			// sort event times
			sampledTimes.sort(Collections.reverseOrder());
			nHiddenEvents += sampledTimes.size();
		}

		currentNode.setHeight(subroot.getHeight());
		currentNode.setID(subroot.getID());
		if (subroot.isLeaf()) {
			currentNode.setNr(subroot.getNr());
		}
		Node firstEvent = putEventsOnBranch(currentNode, sampledTimes);


		if (!subroot.isLeaf()) {
			Node c1 = paintHiddenEvents(subroot.getChild(0), skip);
			Node c2 = paintHiddenEvents(subroot.getChild(1), false);
			currentNode.addChild(c1);
			currentNode.addChild(c2);
		}


		return firstEvent;
	}

	private Node putEventsOnBranch(Node currentNode, List<Double> eventTimes) {
		if (eventTimes.size() == 0) {
			return currentNode;
		}

		Node tmp = new Node();
		Node firstEvent = new Node();
		boolean first = true;
		for (Double eventTime : eventTimes) {
			Node event = new Node();
			if (first) {
				event.setHeight(eventTime);
				firstEvent = event;
				first = false;
			} else {
				tmp.addChild(event);
				event.setHeight(eventTime);
			}
			tmp = event;
		}
		tmp.addChild(currentNode);
		return firstEvent;
	}

	private void updateParameters(int intervalNr) {
		lambda_i = parameterization.getBirthRates()[intervalNr][0];
		mu_i = parameterization.getDeathRates()[intervalNr][0];
		psi_i = parameterization.getSamplingRates()[intervalNr][0];
		t_i = parameterization.getIntervalEndTimes()[intervalNr];
		A_i = A[intervalNr];
		B_i = B[intervalNr];
	}

	private List<Double> sampleTimes(double startTime, double endTime) {
		List<Double> eventTimes = new ArrayList<>();

		double meanN = meanEvents(parameterization.getAge(startTime, finalSampleOffset.getArrayValue()),
				parameterization.getAge(endTime, finalSampleOffset.getArrayValue()));
		if (meanN == 0)
			return eventTimes;

		long n = Randomizer.nextPoisson(meanN);

		for (int s = 0; s < n; s++) {
			eventTimes.add(ftInverse(startTime, endTime, Randomizer.nextDouble()));
		}

		return eventTimes;
	}

	/**
	 * Calculate the mean of the Poisson distribution for the number of hidden
	 * transmission events, that change the host on a branch. It's calculated as
	 * lambda*integral_p_0 at interval i.
	 * 
	 * @param t_0 start of the branch
	 * @param t_1 end of the branch
	 * @return mean of the Poisson distribution for the number of hidden
	 *         transmission events, that change the host on a branch
	 */
	private double meanEvents(double t_0, double t_1) {
		double t0 = t_i - t_0;
		double t1 = t_i - t_1;
		double ans = 0.5 * ((t1 - t0) * (mu_i + psi_i + lambda_i + A_i) + 2.0 * Math
				.log(((-B_i - 1) * Math.exp(A_i * t0) + B_i - 1)
						/ ((-B_i - 1) * Math.exp(A_i * t1) + B_i - 1)));

		return ans;
	}

	private double ft(double t_start, double t_x, double t_end) {
		return meanEvents(parameterization.getAge(t_start, finalSampleOffset.getArrayValue()),
				parameterization.getAge(t_x, finalSampleOffset.getArrayValue())) /
				meanEvents(parameterization.getAge(t_start, finalSampleOffset.getArrayValue()),
						parameterization.getAge(t_end, finalSampleOffset.getArrayValue()));
	}

	private double ftInverse(double t_start, double t_end, double u) {
		double a = 0;
		double b = t_end;
		double b_inf = Double.NEGATIVE_INFINITY;
		double b_sup = Double.POSITIVE_INFINITY;
		for (int j = 0; j < 20; j++) {
			double tmp = (a + b) * 0.5;
			if (ft(t_start, tmp, t_end) <= u) {
				b_inf = tmp;
				b_sup = b;
			}
			if (ft(t_start, tmp, t_end) >= u) {
				b_sup = tmp;
				b_inf = a;
			}
			a = b_inf;
			b = b_sup;
		}
		return (a + b) * 0.5;
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
	
	/**
	 * COPIED from BDMM prime package by T. G. Vaughan on 2022-03-12
	 * 
	 * Apply node numbers to internal nodes below and including subtreeRoot. Numbers
	 * are applied postorder, so parents always have larger numbers than their
	 * children and the root has the hightest number.
	 *
	 * @param subtreeRoot root of subtree
	 * @param nextNumber  next number to be used
	 * @return next number to be used on another part of the tree.
	 */
	private int numberInternalNodesOnSubtree(Node subtreeRoot, int nextNumber) {

		if (subtreeRoot.isLeaf())
			return nextNumber;

		for (Node child : subtreeRoot.getChildren())
			nextNumber = numberInternalNodesOnSubtree(child, nextNumber);

		subtreeRoot.setNr(nextNumber);

		return nextNumber + 1;
	}
	
    /*
     * Loggable implementation
     */

    private long lastRemapSample = -1;

    /**
     * Remap the tree.  Intended to be called by loggers requiring
     * a mapped tree.  Supplying the sample number allows the result
     * of the remapping to be cached and used for other loggers.
     *
     * @param sample sample number at log
     */
    public void remapForLog(long sample) {
        if (!remapOnLogInput.get() || sample == lastRemapSample)
            return;

        map();
        lastRemapSample = sample;
    }

    @Override
    public void init(PrintStream out) {
//		trLog.init(out);
		unmappedTree.init(out);
    }

    @Override
    public void log(long sample, PrintStream out) {
        remapForLog(sample);

//		TransmissionTree tree = (TransmissionTree) getCurrent();
//		out.print("tree STATE_" + sample + " = ");
//        final int[] dummy = new int[1];
//		tree.addOrientationMetadata();
//		String newick = tree.getRoot().toSortedNewick(dummy, true);
//		newick = tree.getRoot().toShortNewick(true);
//        out.print(newick);
//        out.print(";");


		TransmissionTree tree = (TransmissionTree) getCurrent();
		tree.addOrientationMetadata();
		// write out the log tree with meta data
		out.print("tree STATE_" + sample + " = ");
//        tree.getRoot().sort();
//		System.out.println(toNewick(tree.getRoot()));
		out.print(toNewick(tree.getRoot()));
		//out.print(tree.getRoot().toShortNewick(false));
		out.print(";");
    }

    @Override
    public void close(PrintStream out) {
		unmappedTree.close(out);
    }

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
		} //else {
			buf.append(node.getNr() + 1);
//		}

//			if (node.getID() == null) {
//				buf.append(node.getNr() + 1);
//			}
			buf.append("[&");
			buf.append(node.metaDataString);
			buf.append(']');


		buf.append(":");

		double nodeLength;
		nodeLength = node.getLength();

		buf.append(nodeLength);

		return buf.toString();
	}

}
