package tnt.simulator;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
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
	public Input<String> fileNameInput = new Input<>("fileName",
			"Name of file to which gene trees will be written. Only if not using logger.");

	public Input<RealParameter> bottleneckStrengthInput = new Input<RealParameter>("bottleneckStrength",
			"Strength of the bottleneck in scaled time", Validate.REQUIRED);

	public Input<RealParameter> maxBoundInput = new Input<RealParameter>("maxBound",
			"Maximum value of time interval when computing next unobserved transmission time. Default 1000",
			Input.Validate.OPTIONAL);

	public Input<RealParameter> birthRateInput = new Input<RealParameter>("birthRate", "Birth rate",
			Input.Validate.REQUIRED);
	public Input<Function> deathRateInput = new Input<Function>("deathRate", "Death rate", Input.Validate.REQUIRED);
	public Input<RealParameter> samplingRateInput = new Input<RealParameter>("samplingRate",
			"Sampling rate per individual", Input.Validate.REQUIRED);
	public Input<Boolean> completeInput = new Input<>("complete",
			"Is input a complete tree", false);

	public Tree transmissionTree;
    public TraitSet sampleCounts;
	private int transmissionNodeCount;

	private double bottleneckStrength;
	private double maxBound;

	private boolean complete;

	@Override
	public void init(PrintStream out) {
	}

	@Override
	public void close(PrintStream out) {
	}


    @Override
    public void initAndValidate() {
		complete = completeInput.get();

		transmissionTree = transmissionTreeInput.get();
		transmissionNodeCount = transmissionTree.getNodeCount();
		popSizesInput.get().setDimension(transmissionNodeCount);
		bottleneckStrength = bottleneckStrengthInput.get().getArrayValue(0);
		maxBound = maxBoundInput.get() != null ? maxBoundInput.get().getArrayValue(0) : 1000;

        sampleCounts = sampleCountsInput.get();

        assignFromWithoutID(getSimulatedGeneTree());
		// Write simulated network to file if requested
		if (fileNameInput.get() != null) {
			try (PrintStream ps = new PrintStream(fileNameInput.get())) {

				ps.println(toString());

			} catch (FileNotFoundException ex) {
				throw new RuntimeException("Error writing to output file '" + fileNameInput.get() + "'.");
			}
		}
    }

    int getTotalLineageCount(Map<Node, List<Node>> lineages) {
        int count = 0;

        for (List<Node> lineageList : lineages.values())
            count += lineageList.size();

        return count;
    }

    int getTotalSampleCount() {
        int count = 0;

		for (Node speciesNode : transmissionTree.getExternalNodes())
            count += (int)Math.round(sampleCounts.getValue(speciesNode.getID()));

        return count;
    }

    public Tree getSimulatedGeneTree() {

		List<Node> sortedTransmissionTreeNodes = new ArrayList<>(Arrays.asList(transmissionTree.getNodesAsArray()));
		sortedTransmissionTreeNodes.sort(Comparator.comparingDouble(Node::getHeight));

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

			Node minCoal = new Node();
			Node minNonObsTr = new Node();
			double minCoalTime = Double.POSITIVE_INFINITY;
			double minNonObsTrTime = Double.POSITIVE_INFINITY;

			for (Node transmissionNode : activeLineages.keySet()) {
				int k = activeLineages.get(transmissionNode).size();
				if (k > 1) {
					double timeToNextCoal = Randomizer
							.nextExponential(0.5 * k * (k - 1) * 1 / popSize.get(transmissionNode));
					if (timeToNextCoal < minCoalTime) {
						minCoalTime = timeToNextCoal;
						minCoal = transmissionNode;
					}

					if (!complete) {
					double timeToNextNonObsTr = soveForTime(t,
							birthRateInput.get().getArrayValue(0), deathRateInput.get().getArrayValue(0),
							samplingRateInput.get().getArrayValue(0));
					if (timeToNextNonObsTr < minNonObsTrTime) {
						minNonObsTrTime = timeToNextNonObsTr;
						minNonObsTr = transmissionNode;

					}
					}
				}
			}

			
			double dt = Math.min(minCoalTime, minNonObsTrTime);

			if (!sortedTransmissionTreeNodes.isEmpty() && t + dt > sortedTransmissionTreeNodes.get(0).getHeight()) {
				Node transmissionNode = sortedTransmissionTreeNodes.get(0);
				t = transmissionNode.getHeight();

				activeLineages.put(transmissionNode, new ArrayList<>());

				if (transmissionNode.isLeaf()) {
				    // Sample

					int count = (int) Math.round(sampleCounts.getValue(transmissionNode.getID()));

                    for (int i=0; i<count; i++) {
						Node geneTreeSampleNode = new Node(transmissionNode.getID() + (i + 1));
                        geneTreeSampleNode.setNr(nextLeafNodeNr++);
						geneTreeSampleNode.setHeight(transmissionNode.getHeight());
						activeLineages.get(transmissionNode).add(geneTreeSampleNode);
                    }

                } else {
					// Observed transmission

					Node recipient = transmissionNode.getChildren().get(1);
					List<Node> recipientLineages = activeLineages.get(recipient);
					double recipientNe = popSize.get(recipient);

					// From Tim: Can this be extracted to a method?  More-or-less duplicated below.
					if (bottleneckStrength > 0 && recipientLineages.size() > 1 && !recipient.isDirectAncestor()) {
						double duplicateTime = t;
						double stopTime = t + bottleneckStrength;
						while (duplicateTime < stopTime) {
							int nrLineages = recipientLineages.size();
							double deltaT = Randomizer
									.nextExponential(nrLineages * (nrLineages - 1) * 0.5 * 1 / recipientNe);
							if (duplicateTime + deltaT < stopTime) {
								int k = recipientLineages.size();
								Node node1 = recipientLineages.get(Randomizer.nextInt(k));
								Node node2;
								do {
									node2 = recipientLineages.get(Randomizer.nextInt(k));
								} while (node2 == node1);

								Node parent = new Node(String.valueOf(nextIntNodeNr));
								parent.setNr(nextIntNodeNr++);
								parent.setHeight(t);
								parent.addChild(node1);
								parent.addChild(node2);
								recipientLineages.remove(node1);
								recipientLineages.remove(node2);
								recipientLineages.add(parent);
							} else {
								break;
							}
							duplicateTime += deltaT;
						}
						activeLineages.put(recipient, recipientLineages);
					}

							
							
					for (Node speciesChild : transmissionNode.getChildren()) {

						activeLineages.get(transmissionNode).addAll(activeLineages.get(speciesChild));
						activeLineages.get(speciesChild).clear();
					}
                }

				sortedTransmissionTreeNodes.remove(0);

			} else if (dt != Double.POSITIVE_INFINITY) {
			    // Coalescence

                t += dt;

				if (minCoalTime == dt) {
					List<Node> lineageList = activeLineages.get(minCoal);
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
					activeLineages.put(minCoal, lineageList);
					
				} else if (minNonObsTrTime == dt) {
					// Unobserved transmission

					List<Node> lineageList = activeLineages.get(minNonObsTr);
					double Ne = popSize.get(minNonObsTr);

					// Duplicate of above code
					if (bottleneckStrength > 0 && lineageList.size() > 1) {
						double duplicateTime = t;
						double stopTime = t + bottleneckStrength;
						while (duplicateTime < stopTime) {
							int nrLineages = lineageList.size();
							double deltaT = Randomizer
									.nextExponential(nrLineages * (nrLineages - 1) * 0.5 * 1 / Ne);
							if (duplicateTime + deltaT < stopTime) {
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
							} else {
								break;
							}
							duplicateTime += deltaT;
						}
						activeLineages.put(minNonObsTr, lineageList);
					}


				}
            }
        }

        // Return tree with remaining lineage as root
		return new Tree(activeLineages.get(transmissionTree.getRoot()).get(0));
    }

	// Solve an inverse transform problem for time t of next occurence
	private double soveForTime(double currentTime, double lambda, double mu, double psi) {

		double u = Randomizer.nextDouble();
		double c_1 = Math.abs(Math.sqrt( Math.pow((lambda - mu - psi), 2) + 4*lambda*psi ));
		double c_2 = - (lambda - mu - psi)/c_1;
		
		NewtonRaphsonSolver solver = new NewtonRaphsonSolver(1E-10);

		UnivariateDifferentiableFunction f = new UnivariateDifferentiableFunction() {

			private double function(double t) {
				return ((-2 * Math.log(Math.abs((c_2 - 1) * Math.exp(-c_1 * t) - c_2 - 1))
						+ (-c_1 * t) + (mu + psi + lambda) * t) / 2);

			}

			private DerivativeStructure ds(DerivativeStructure t) {
				return ((t.multiply(-c_1)).exp().multiply(c_2 - 1).subtract(c_2 + 1).abs().log().multiply(-2)
						.add(t.multiply(-c_1).exp().log()).add(t.multiply(lambda + mu + psi))).divide(2);
			}

			// solve fot t: int_{currentTime}^{t} f(x) = -ln(u)
			@Override
			public double value(double t) {
				return this.function(t) + Math.log(u) - this.function(currentTime);
			}

			@Override
			public DerivativeStructure value(DerivativeStructure t) throws DimensionMismatchException {

				return this.ds(t).add(Math.log(1 - u)).subtract(this.ds(new DerivativeStructure(1, 1, currentTime)));

			}
		};

		return solver.solve(10000, f, currentTime, maxBound);

	}

}
