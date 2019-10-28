package tnt.simulator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.analysis.differentiation.UnivariateDifferentiableFunction;
import org.apache.commons.math3.analysis.solvers.NewtonRaphsonSolver;
import org.apache.commons.math3.exception.DimensionMismatchException;

import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

/**
 * @author Ugne Stolz. Created on 22 Oct 2019 Adapted from Tim Vaughan's
 *         SimulatedGenehTree for gene trees conditioned on species tree in
 *         StarBeast2: https://doi.org/10.1093/molbev/msx126
 */

public class SimulatedGeneTree extends Tree {

	public Input<Tree> transmissionTreeInput = new Input<>("transmissionTreeInput",
			"Fully labeled transmission tree on which to simulate gene trees",
            Input.Validate.REQUIRED);

    public Input<TraitSet> sampleCountsInput = new Input<>(
            "sampleCounts",
            "TraitSet defining number of  samples per node in species tree.",
            Input.Validate.REQUIRED);
// For now we assume constant population sizes along each branch
	public Input<RealParameter> popSizesInput = new Input<RealParameter>("populationSizes",
			"Constant per-branch effective population sizes.", Validate.REQUIRED);

	public Input<RealParameter> bottleneckStrength = new Input<RealParameter>("bottleneckStrength",
			"Strength of the bottleneck in scaled time", Validate.REQUIRED);

	public Input<RealParameter> birthRateInput = new Input<RealParameter>("birthRate", "Birth rate",
			Input.Validate.REQUIRED);
	public Input<Function> deathRateInput = new Input<Function>("deathRate", "Death rate", Input.Validate.REQUIRED);
	public Input<RealParameter> samplingRateInput = new Input<RealParameter>("samplingRate",
			"Sampling rate per individual", Input.Validate.REQUIRED);
	// Add rho this to P_0
	public Input<RealParameter> rhoProbability = new Input<RealParameter>("rho",
			"Probability of an individual to be sampled at present", (RealParameter) null);

	public Tree transmissionTree;
    public TraitSet sampleCounts;
	private int transmissionNodeCount;

    public SimulatedGeneTree() { }

    @Override
    public void initAndValidate() {
		transmissionTree = transmissionTreeInput.get();
		transmissionNodeCount = transmissionTree.getNodeCount();
		popSizesInput.get().setDimension(transmissionNodeCount);

        sampleCounts = sampleCountsInput.get();

        assignFromWithoutID(getSimulatedGeneTree());
    }

    int getTotalLineageCount(Map<Node, List<Node>> lineages) {
        int count = 0;

        for (List<Node> lineageList : lineages.values())
            count += lineageList.size();

        return count;
    }

//    TODO what about the sampling through time?
    int getTotalSampleCount() {
        int count = 0;

		for (Node speciesNode : transmissionTree.getExternalNodes())
            count += (int)Math.round(sampleCounts.getValue(speciesNode.getID()));

        return count;
    }

    public Tree getSimulatedGeneTree() {

		List<Node> sortedTransmissionTreeNodes = new ArrayList<>(Arrays.asList(transmissionTree.getNodesAsArray()));

		sortedTransmissionTreeNodes.sort((o1, o2) -> {
            if (o1.getHeight() < o2.getHeight())
                return -1;
            if (o1.getHeight() > o2.getHeight())
                return 1;

            return 0;
        });

        // Perform simulation

        Map<Node,List<Node>> activeLineages = new HashMap<>();
		Map<Node, Double> popSize = new HashMap<>();
        int nextLeafNodeNr = 0;
        int nextIntNodeNr = getTotalSampleCount();
        double t = 0.0;

		int ind = 0;
		for (Node transmissionNode : sortedTransmissionTreeNodes) {
			popSize.put(transmissionNode, popSizesInput.get().getArrayValue(ind));
			ind++;
		}

		while (getTotalLineageCount(activeLineages) > 1 || !sortedTransmissionTreeNodes.isEmpty()) {

            // Compute propensity

            double totalPropensity = 0;
            Map<Node, Double> propensities = new HashMap<>();
            Map<Node, Double> nonObsTransmissionPropensities = new HashMap<>();
			for (Node transmissionNode : activeLineages.keySet()) {
				int k = activeLineages.get(transmissionNode).size();
				double thisProp = 0.5 * k * (k - 1) * 1 / popSize.get(transmissionNode);
				propensities.put(transmissionNode, thisProp);
                totalPropensity += thisProp;
                
				double nonObservedTransmission = getInverseIntensity(t, transmissionNode.getParent().getHeight(),
						birthRateInput.get().getArrayValue(0), deathRateInput.get().getArrayValue(0),
						samplingRateInput.get().getArrayValue(0));
				nonObsTransmissionPropensities.put(transmissionNode, nonObservedTransmission);
            }

			double dt_coal = Randomizer.nextExponential(totalPropensity);
			double dt_nonObservedCoal = Collections.min(nonObsTransmissionPropensities.values());
			
			double dt = Math.min(dt_coal, dt_nonObservedCoal);

			if (!sortedTransmissionTreeNodes.isEmpty() && t + dt > sortedTransmissionTreeNodes.get(0).getHeight()) {
				Node transmissionNode = sortedTransmissionTreeNodes.get(0);
				t = transmissionNode.getHeight();

				activeLineages.put(transmissionNode, new ArrayList<>());

				if (transmissionNode.isLeaf()) {
					int count = (int) Math.round(sampleCounts.getValue(transmissionNode.getID()));

                    for (int i=0; i<count; i++) {
						Node geneTreeSampleNode = new Node(transmissionNode.getID() + String.valueOf(i + 1));
                        geneTreeSampleNode.setNr(nextLeafNodeNr++);
						geneTreeSampleNode.setHeight(transmissionNode.getHeight());
						activeLineages.get(transmissionNode).add(geneTreeSampleNode);
                    }

                } else {
					Node child_1 = transmissionNode.getChildren().get(0);
					Node child_2 = transmissionNode.getChildren().get(1);

					String donorID = transmissionNode.getID();
					if (donorID != child_1.getID() || donorID != child_2.getID()) {
						// TODO exception
						System.out.println("DEBUG: Wrongly labeled transmission tree");
						System.exit(0);
					}
					String recipientID = child_1.getID() != donorID ? child_1.getID() : child_2.getID();
							
							
//					for (Node speciesChild : transmissionNode.getChildren()) {
//						String donorID = transmissionNode.getID();
//						String recipientID = transmissionNode.getID() !=
//						
//						activeLineages.get(transmissionNode).addAll(activeLineages.get(speciesChild));
//                        activeLineages.get(speciesChild).clear();
//                    }
                }

				sortedTransmissionTreeNodes.remove(0);

            } else {
				// TODO here need to include two variants
                t += dt;

                // Coalesce a random pair of lineages

                double u = Randomizer.nextDouble()*totalPropensity;
				for (Node transmissionNode : propensities.keySet()) {
					u -= propensities.get(transmissionNode);
                    if (u < 0) {
						List<Node> lineageList = activeLineages.get(transmissionNode);
                        int k = lineageList.size();

                        Node node1 = lineageList.get(Randomizer.nextInt(k));
                        Node node2;
                        do {
                            node2 = lineageList.get(Randomizer.nextInt(k));
                        } while (node2 == node1);

                        Node parent = new Node(String.valueOf(nextIntNodeNr));
                        parent.setNr(nextIntNodeNr++);
                        parent.setHeight(t);
                        parent.addChild(node1);
                        parent.addChild(node2);
                        lineageList.remove(node1);
                        lineageList.remove(node2);
                        lineageList.add(parent);

                        break;
                    }
                }
            }
        }

        // Return tree with remaining lineage as root
		return new Tree(activeLineages.get(transmissionTree.getRoot()).get(0));
    }

	private double getInverseIntensity(double currentTime, double parentTime, double lambda, double mu, double psi) {
		double u = Randomizer.nextDouble();
		double c_1 = Math.abs(Math.sqrt( Math.pow((lambda - mu - psi), 2) + 4*lambda*psi ));
		double c_2 = - (lambda - mu - psi)/c_1;
		
		NewtonRaphsonSolver solver = new NewtonRaphsonSolver(1E-10);

		UnivariateDifferentiableFunction f = new UnivariateDifferentiableFunction() {

			private double function(double t) {
				return ((-2 * Math.log(Math.abs((c_2 - 1) * Math.exp(-c_1 * t) - c_2 - 1))
						+ Math.log(Math.exp(-c_1 * t)) + (mu + psi + lambda) * t) / 2);

			}

			private DerivativeStructure ds(DerivativeStructure t) {
				return ((t.multiply(-c_1)).exp().multiply(c_2 - 1).subtract(c_2 + 1).abs().log().multiply(-2)
						.add(t.multiply(-c_1).exp().log()).add(t.multiply(lambda + mu + psi))).divide(2);
			}

			@Override
			public double value(double t) {
				System.out.println(this.function(currentTime));
				return this.function(t) + Math.log(1 - u) - this.function(currentTime);
			}

			@Override
			public DerivativeStructure value(DerivativeStructure t) throws DimensionMismatchException {

				return this.ds(t).add(Math.log(1 - u)).subtract(this.ds(new DerivativeStructure(1, 1, currentTime)));

			}
		};

		return solver.solve(1000, f, currentTime, parentTime);

	}

}
