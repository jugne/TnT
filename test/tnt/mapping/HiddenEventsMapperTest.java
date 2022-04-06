package tnt.mapping;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.Test;

import bdmmprime.parameterization.CanonicalParameterization;
import bdmmprime.parameterization.Parameterization;
import bdmmprime.parameterization.SkylineVectorParameter;
import bdmmprime.parameterization.TimedParameter;
import bdmmprime.parameterization.TypeSet;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import org.junit.Assert;
import tnt.simulator.SABDSimulator;

import org.apache.commons.math3.stat.inference.TestUtils;

public class HiddenEventsMapperTest {

	double lambda = 3.57383023240944; // birth/transmission rate
	double mu = 1.12242310452165; // death rate
	double psi = 3.41495738932089; // sampling rate
	double r = 0.5; // removal upon sampling probability
	double rho = 1.0; // sampling at present prob
	double origin = 1.0; // origin of the transmission process
	int minSamples = 1; // minimum samples for tree to be accepted after simulation
	int maxSamples = 2000; // maximum samples for tree to be accepted after simulation
	int nSims = 100000; // number of simulations
	int tmpHiddenSim;

	@Test
	public void mapperTest() throws Exception {

		long seed = Randomizer.nextLong();
		System.out.println(seed);
		Randomizer.setSeed(seed);

		int[] hiddenEventsSim = new int[nSims];
		int[] hiddenEventsMap = new int[nSims];

		SABDSimulator simulator = new SABDSimulator(lambda,
				mu, psi, r, rho, origin);

		for (int i = 0; i < nSims; i++) {

			Randomizer.setSeed(Randomizer.nextLong());

			int result;
			do {
				result = simulator.simulate();
				if (result < 0 || simulator.sampledNodeNumber < minSamples
						|| simulator.sampledNodeNumber > maxSamples) {
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

			tmpHiddenSim = 0;
			updateHiddenEventsCounter(simulator.fullRoot, sampledTreeLeafNodeIds);
			hiddenEventsSim[i] = tmpHiddenSim;

			MappedTree mappedTree = new MappedTree();
			Parameterization parameterization = new CanonicalParameterization();
			parameterization.initByName(
					"origin", new RealParameter(Double.toString(origin)),
					"typeSet", new TypeSet(1),
					"birthRate", new SkylineVectorParameter(
							null,
							new RealParameter(Double.toString(lambda))),
					"deathRate", new SkylineVectorParameter(
							null,
							new RealParameter(Double.toString(mu))),
					"samplingRate", new SkylineVectorParameter(
							null,
							new RealParameter(Double.toString(psi))),
					"removalProb", new SkylineVectorParameter(
							null,
							new RealParameter(Double.toString(r))),
					"rhoSampling", new TimedParameter(
							new RealParameter(Double.toString(origin)),
							new RealParameter(Double.toString(rho))));

			Tree t = new Tree(simulator.sampledRoot.toNewick());

			Double offset = origin + smallestHeightSampled;

			mappedTree.initByName(
					"tree", t,
					"parameterization", parameterization,
					"finalSampleOffset", new RealParameter(offset.toString()));

			IntegerParameter h = (IntegerParameter) mappedTree.hiddenEventsCounterInput.get();
			hiddenEventsMap[i] = h.getValue();

		}

		int maxSim = Arrays.stream(hiddenEventsSim).max().getAsInt();
		int maxMap = Arrays.stream(hiddenEventsMap).max().getAsInt();

		int maxMax = Math.max(maxSim, maxMap);

		int[] histMap = calcHistogram(hiddenEventsMap, 0.5, maxMax + 0.5, maxMax);
		int[] histSim = calcHistogram(hiddenEventsSim, 0.5, maxMax + 0.5, maxMax);

		List<Integer> hiddenEventsSimList = new ArrayList<Integer>();
		List<Integer> hiddenEventsMapList = new ArrayList<Integer>();

		for (int j = 0; j < histMap.length; j++) {
			if (histMap[j] != 0 || histSim[j] != 0) {
				hiddenEventsSimList.add(histSim[j]);
				hiddenEventsMapList.add(histMap[j]);
			}
		}

		double[] hiddenEventsSimDouble = hiddenEventsSimList.stream()
				.mapToDouble(i -> i).toArray();
		long[] hiddenEventsMapLong = hiddenEventsMapList.stream().mapToLong(i -> i).toArray();

		// https://commons.apache.org/proper/commons-math/userguide/stat.html
//		System.out.println(TestUtils.chiSquareTest(hiddenEventsSimDouble,
//				hiddenEventsMapLong));
		Assert.assertTrue(TestUtils.chiSquareTest(hiddenEventsSimDouble,
				hiddenEventsMapLong, 0.01));
		;

	}

	public static int[] calcHistogram(int[] data, double min, double max, int numBins) {
		final int[] result = new int[numBins];
		final double binSize = (max - min) / numBins;

		for (int d : data) {
			int bin = (int) ((d - min) / binSize);
			if (bin < 0) {
				/* this data is smaller than min */ } else if (bin >= numBins) {
				/* this data point is bigger than max */ } else {
				result[bin] += 1;
			}
		}
		return result;
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
//					hidden = false;
					}
				}
				if (hidden) {
					tmpHiddenSim += 1;
				}

			}
			updateHiddenEventsCounter(fullTreeSubroot.getRight(), sampledTreeLeafNodeIds);
			updateHiddenEventsCounter(fullTreeSubroot.getLeft(), sampledTreeLeafNodeIds);
		}
	}

}
