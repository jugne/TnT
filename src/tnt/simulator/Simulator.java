package tnt.simulator;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import beast.evolution.tree.Node;
import beast.util.Randomizer;

/**
 * @author originally Alexandra Gavryushkina, adapted by Ugne Stolz
 *         <ugne.stolz@protonmail.com>
 * @date 12 Mar 2022
 */

public class Simulator {

	double[] parameters = new double[9];

    public Simulator(String[] parametersStr) {
		for (int i = 0; i < 9; i++) {
            parameters[i] = Double.parseDouble(parametersStr[i]); 
        }
    }

    public static void main(String[] args) throws Exception {

		if (args.length != 9) {
			System.out.println(
					"There have to be 9 arguments for parameters: d (lambda), nu (mu), s (psi), r, rho, t_origin, seed, minSamples, maxSamples");
        } else {
            Simulator simulator = new Simulator(args);
            simulator.simulateForTotalEvidence();
        }

    }

    private void simulateForTotalEvidence() throws Exception {

        PrintStream writer = null;
//        PrintStream treeWriter = null;

        int treeCount = 1;

        int index=0;
        int count=0;


        try {
            writer = new PrintStream(new File("trees_and_pars.txt"));
//            treeWriter = new PrintStream(new File("trees.txt"));


            for (int i=0; i< treeCount; i++) {
                //int meanLeafCount = 0;
                //int low=0;
                //int high=0;

                //parameters = {0.04, 0.6, 0.2, 0.0, 0.8};
//				SABDSimulator simulator = new SABDSimulator(1.5217, 0.5, 0.7575, 0.5, 0.0, true, 3.0);
				Randomizer.setSeed((long) parameters[6]);
				SABDSimulator simulator = new SABDSimulator(parameters[0], parameters[1], parameters[2], parameters[3],
						parameters[4], parameters[5]);
                int result;
                do {
                    result = simulator.simulate();
					if (result < 0 || simulator.sampledNodeNumber < parameters[7]
							|| simulator.sampledNodeNumber > parameters[8]) {
						result = -1;
						count++;
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

				writer.println("full tree");
				writer.println(simulator.fullRoot.toNewick());
				writer.println("sampled tree");
                writer.println(simulator.sampledRoot.toNewick());
                //treeWriter.println(simulator.root.toShortNewick(false) + ";");

				writer.println("full tree traits");
				simulator.printTraitsWithRhoSamplingTime(simulator.fullRoot, writer, 0.0);
				writer.println("sampled tree traits");
                simulator.printTraitsWithRhoSamplingTime(simulator.sampledRoot, writer, (simulator.rhoSamplingTime + smallestHeightSampled));

				writer.println("same patient samples");
				simulator.printSamepatientSamples(simulator.sampledRoot, writer);

				writer.println("origin");
				writer.println(simulator.rhoSamplingTime);


				writer.println("lambda");
				writer.println(parameters[0]);
				writer.println("mu");
				writer.println(parameters[1]);
				writer.println("psi");
				writer.println(parameters[2]);
				writer.println("r");
				writer.println(parameters[3]);
				writer.println("rho");
				writer.println(parameters[4]);
				writer.println("number of hosts");
				writer.println(simulator.hosts);
				writer.println("transmission times");
				simulator.printTransmissionTimes(simulator.fullRoot, writer);
				writer.println("sampled transmission times");
				simulator.printTransmissionTimes(simulator.fullRoot, writer);
				writer.println("SA count");
				writer.println(simulator.countSA(simulator.sampledRoot));
				writer.println("past sample count");
				writer.println(simulator.sampledNodeNumber - simulator.rhoSampledNodeNumber);
				writer.println("rho sample count");
				writer.println(simulator.rhoSampledNodeNumber);
				writer.println("total sample count");
				writer.println(simulator.sampledNodeNumber);
				writer.println("offset");
				writer.println(simulator.rhoSamplingTime + smallestHeightSampled);
//				writer.println(simulator.rhoSamplingTime + smallestHeightSampled);
				writer.println("number rejected simulations:" + count);

            }
            //System.out.print(meanLeafCount/(treeCount+count));
            //System.out.println();
            System.out.println("Number of trees rejected: " + count);


        } catch (IOException e) {
			System.out.println(e.getMessage());
        }
		finally {
			if (writer != null) {
				writer.close();
			}
//            if (treeWriter != null) {
//                treeWriter.close();
//            }
		}

    }
}
