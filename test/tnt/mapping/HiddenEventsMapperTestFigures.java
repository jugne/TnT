package tnt.mapping;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import bdmmprime.parameterization.CanonicalParameterization;
import bdmmprime.parameterization.Parameterization;
import bdmmprime.parameterization.SkylineMatrixParameter;
import bdmmprime.parameterization.SkylineVectorParameter;
import bdmmprime.parameterization.TimedParameter;
import bdmmprime.parameterization.TypeSet;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import tnt.mapping.MappedTree;
import tnt.simulator.SABDSimulator;

/**
 * @author original simulator Alexandra Gavryushkina, adapted by for hidden
 *         events counting Ugne Stolz <ugne.stolz@protonmail.com>
 * @date 12 Mar 2022
 */

public class HiddenEventsMapperTestFigures {

	int hiddenEventsCounter;
	double rhoSamplingTime;

	PrintStream writerHiddenEventsTime;

	String outputName;

	double[] parameters = new double[10];

	public HiddenEventsMapperTestFigures(String[] parametersStr) {
		for (int i = 0; i < 10; i++) {
			parameters[i] = Double.parseDouble(parametersStr[i]);
		}
		outputName = parametersStr[10];
	}

	public static void main(String[] args) throws Exception {

		if (args.length != 11) {
			System.out.println(
					"There have to be 10 arguments for parameters: lambda, mu, psi, r, rho, t_origin, seed, minSamples, maxSamples, nSims, outputName");
		} else {

			HiddenEventsMapperTestFigures simulator = new HiddenEventsMapperTestFigures(args);
			simulator.simulateForTotalEvidence();
		}

	}

	private void simulateForTotalEvidence() throws Exception {

		PrintStream writerFullTree = null;
		PrintStream writerSampledTrees = null;
		PrintStream writerMappedTrees = null;
		PrintStream writerHiddenEvents = null;
		PrintStream writerHiddenEventsMap = null;
		writerHiddenEventsTime = null;

		try {

			writerFullTree = new PrintStream(new File(outputName + ".fullTree.trees"));
			writerSampledTrees = new PrintStream(new File(outputName + ".sampledTree.trees"));
			writerMappedTrees = new PrintStream(new File(outputName + ".mappedTree.trees"));
			writerHiddenEventsMap = new PrintStream(new File(outputName + ".hiddenEventsMap.log"));
			writerHiddenEventsTime = new PrintStream(new File(outputName + ".hiddenEventsTimes.log"));
			writerHiddenEventsMap.println("sample" + "\t" + "hiddenEventsSim" + "\t" + "hiddenEventsMap");
			writerHiddenEventsTime.println("sample" + "\t" + "hiddenEventsTimes");

			if ((long) parameters[6] < 0) {
				Randomizer.setSeed(Randomizer.nextLong());
			} else
				Randomizer.setSeed((long) parameters[6]);
			rhoSamplingTime = parameters[4];

			SABDSimulator simulator = new SABDSimulator(parameters[0],
					parameters[1], parameters[2], parameters[3],
					parameters[4], parameters[5]);


			for (int i = 0; i < parameters[9]; i++) {

				Randomizer.setSeed(Randomizer.nextLong());
				hiddenEventsCounter = 0;
				int result;
				do {
					result = simulator.simulate();
					if (result < 0 || simulator.sampledNodeNumber < parameters[7]
							|| simulator.sampledNodeNumber > parameters[8]) {
						result = -1;
					}
				} while (result < 0);

				double smallestHeightSampled = Double.POSITIVE_INFINITY;
				for (Node n : simulator.sampledRoot.getAllLeafNodes()) {
					if (n.getHeight() < smallestHeightSampled)
						smallestHeightSampled = n.getHeight();
				}
				double smallestHeightFull = Double.POSITIVE_INFINITY;
				for (Node n : simulator.fullRoot.getAllLeafNodes()) {
					if (n.getHeight() < smallestHeightFull)
						smallestHeightFull = n.getHeight();
				}

				List<Node> sampledTreeLeafNodes = simulator.sampledRoot.getAllLeafNodes();
				List<String> sampledTreeLeafNodeIds = new ArrayList<String>();

				for (Node k : sampledTreeLeafNodes) {
					sampledTreeLeafNodeIds.add(k.getID());

				}
				writerHiddenEventsTime.print(i);
				updateHiddenEventsCounter(simulator.fullRoot, sampledTreeLeafNodeIds);

				writerFullTree.println(simulator.fullRoot.toNewick() + ";");

				writerSampledTrees.println(simulator.sampledRoot.toNewick() + ";");

				MappedTree mappedTree = new MappedTree();
				Parameterization parameterization = new CanonicalParameterization();
				parameterization.initByName(
						"origin", new RealParameter(Double.toString(parameters[5])),
						"typeSet", new TypeSet(1),
						"birthRate", new SkylineVectorParameter(
								null,
								new RealParameter(Double.toString(parameters[0]))),
						"deathRate", new SkylineVectorParameter(
								null,
								new RealParameter(Double.toString(parameters[1]))),
						"samplingRate", new SkylineVectorParameter(
								null,
								new RealParameter(Double.toString(parameters[2]))),
						"removalProb", new SkylineVectorParameter(
								null,
								new RealParameter(Double.toString(parameters[3]))),
						"rhoSampling", new TimedParameter(
								new RealParameter(Double.toString(parameters[5])),
								new RealParameter(Double.toString(parameters[4]))));

				Tree t = new Tree(simulator.sampledRoot.toNewick());

				Double offset = parameters[5] + smallestHeightSampled;

				mappedTree.initByName(
						"tree", t,
						"parameterization", parameterization,
						"finalSampleOffset", new RealParameter(offset.toString()));

				IntegerParameter h = (IntegerParameter) mappedTree.hiddenEventsCounterInput.get();

				writerHiddenEventsMap.println(i + "\t" + hiddenEventsCounter + "\t" + h.getValue());
				writerMappedTrees.println(mappedTree.getRoot().toNewick() + ";");
			}
		} catch (IOException e) {
			System.out.println(e.getMessage());
		} finally {
			if (writerFullTree != null) {
				writerFullTree.close();
			}
			if (writerSampledTrees != null) {
				writerSampledTrees.close();
			}
			if (writerHiddenEvents != null) {
				writerHiddenEvents.close();
			}
			if (writerHiddenEventsMap != null) {
				writerHiddenEventsMap.close();
			}
		}
	}

	private void updateHiddenEventsCounter(Node fullTreeSubroot, List<String> sampledTreeLeafNodeIds) {
		if (!fullTreeSubroot.isLeaf()) {
			boolean hidden = false;
			if (!fullTreeSubroot.isFake()) {
				List<Node> fullTreeLeavesRight = new ArrayList<>();
				fullTreeSubroot.getRight().getAllLeafNodes(fullTreeLeavesRight);
				List<Node> fullTreeLeavesLeft = new ArrayList<>();
				fullTreeSubroot.getLeft().getAllLeafNodes(fullTreeLeavesLeft);

				for (Node n : fullTreeLeavesRight) {
					if (sampledTreeLeafNodeIds.contains(n.getID())) {
						hidden = true;
						break;
					}
				}

				if (hidden) {
					for (Node n : fullTreeLeavesLeft) {
						if (sampledTreeLeafNodeIds.contains(n.getID())) {
							hidden = false;
							break;
						}
//						hidden = false;
					}
				}
				if (hidden) {
					hiddenEventsCounter += 1;
					writerHiddenEventsTime.print("\t" + (rhoSamplingTime
							+ fullTreeSubroot.getHeight()));
				}

			}
			updateHiddenEventsCounter(fullTreeSubroot.getRight(), sampledTreeLeafNodeIds);
			updateHiddenEventsCounter(fullTreeSubroot.getLeft(), sampledTreeLeafNodeIds);
		}
	}

}