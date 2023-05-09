package tnt.tntAnnotator;

import beast.app.treeannotator.TreeAnnotator;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import tnt.transmissionTree.TransmissionTree;
import tnt.util.Tools;

import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.List;
import java.util.*;

//TODO make burnin percentage optional

public class IndirectTransmissionAnalyser2 extends TreeAnnotator {

	private static class TntAnalyserOptions extends TreeAnnotator {
		File inTransmissionTreeFile;
		File outFile;

		@Override
		public String toString() {
			return "Active options:\n" +
					"Input transmission tree file: " + inTransmissionTreeFile + "\n" +
					"Output file: " + outFile + "\n";
		}
	}

	public class SeparatingBifurcations {
		public List<String> taxaList;
		public int n_transmissions;

		public SeparatingBifurcations(List<String> taxaList, int n_transmissions) {
			this.taxaList = taxaList;
			this.n_transmissions = n_transmissions;
		}
	}

	public IndirectTransmissionAnalyser2(TntAnalyserOptions options) throws IOException {
 	
    	
        // Display options:
        System.out.println(options + "\n");

        // Initialise reader

		TreeSet trTreeSet = new FastTreeSet(options.inTransmissionTreeFile.toString(), 10);
		// new MemoryFriendlyTreeSet(options.inTransmissionTreeFile.toString(), 0);//
		//		Set<String> leafNodes = new HashSet<String>(options.leafIds);

		PrintStream ps = new PrintStream(options.outFile);
		boolean first = true;
		List<String> hostsList = null;


		trTreeSet.reset();
		while (trTreeSet.hasNext()) {
			Tree tree = trTreeSet.next();
			Set<String> hosts = new HashSet<String>();
			if (first) {
				hosts.add("0");
				for (int i=0; i<tree.getLeafNodeCount(); i++) {
					hosts.add(tree.getNode(i).getID().split("_")[0]);
				}
				hostsList = new ArrayList<>(hosts);
				for (int i = 0; i < hostsList.size(); i++) {
					for (int j = 0; j < hostsList.size(); j++) {
						if (i == hostsList.size() - 1 && j == hostsList.size() - 1)
							ps.print(hostsList.get(i) + "_" + hostsList.get(j));
						else
							ps.print(hostsList.get(i) + "_" + hostsList.get(j) + "\t");
					}
				}
				ps.print("\n");
			}
			first = false;
			Double[] transmissions = new Double[hostsList.size() * hostsList.size()];
			Arrays.fill(transmissions, 0.0);
			String n = tree.getRoot().toNewick();
			for (Node leaf : tree.getExternalNodes()){
				assignIds(leaf);
			}
			n = tree.getRoot().toNewick();
			for (Node leaf : tree.getExternalNodes()){
				fillTransmissions(leaf, transmissions, hostsList);
			}
//			BitSet transmissions = new BitSet(hostsList.size() * hostsList.size());
//			HashMap<String, List<Double>> transmissionTimes = new HashMap<String, List<Double>>();
//			fillTransmissions(tree.getRoot(), transmissions, transmissionTimes, hostsList);


			for (int i = 0; i < hostsList.size() * hostsList.size(); i++) {
				if (i == hostsList.size() * hostsList.size() - 1) {
					ps.print(transmissions[i]);
				} else {
					ps.print(transmissions[i] + "\t");
				}
			}
			ps.print("\n");
        }
		System.out.println("\nDone!");
 

        }

		private void assignIds(Node leaf){
			if (leaf.isRoot())
				return;
			String id = leaf.getID().split("_")[0];
			Node parent = leaf.getParent();

			String parentId = "";
			if (parent.isFake()){
				parentId = Tools.equalHeightWithPrecision(parent,leaf) ? id : getOtherChild(leaf).getID().split("_")[0];
				parent.setID(parentId+"_");
				assignIds(parent);
			} else if (parent.getChildCount()==1){
				return;
			} else if (leaf.metaDataString.contains("orientation=donor")){
				parent.setID(id+"_");
				assignIds(parent);
			}
		}

	private void fillTransmissions(Node leaf, Double[] transmissions,
								   List<String> hostsList){
		String recipientID = leaf.getID().split("_")[0];
		Node parent = leaf.getParent();
		while(!parent.isRoot() &&
				(parent.getChildCount() == 1 ||
						Objects.equals(parent.getID().split("_")[0], leaf.getID().split("_")[0]) ||
						!parent.getID().contains("_"))){
			parent = parent.getParent();
		}
		String donorID = parent.getID().split("_")[0];
		int idx = hostsList.indexOf(donorID);
		if (idx<0){
			idx = hostsList.indexOf("0");
			idx = idx * hostsList.size() + hostsList.indexOf(recipientID);
			transmissions[idx] = 1.0;
			return;
		}
		idx = idx * hostsList.size() + hostsList.indexOf(recipientID);
		transmissions[idx] = 1.0;


	}

	private Node getOtherChild(Node child){
			Node parent = child.getParent();
			return parent.getChild(0).getNr()==child.getNr() ? parent.getChild(1) : parent.getChild(0);
	}



	/**
	 * Prepare JFrame to which tntOrientator output streams will be directed.
	 */
	private static void setupGUIOutput() {

		JFrame frame = new JFrame();
		frame.setTitle("Transmission Tree Orientator");
		frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);

		JTextArea textArea = new JTextArea(25, 80);
		textArea.setFont(new Font("monospaced", Font.PLAIN, 12));
		textArea.setEditable(false);
		frame.getContentPane().add(new JScrollPane(textArea), BorderLayout.CENTER);

		JButton closeButton = new JButton("Close");
		closeButton.addActionListener(e -> System.exit(0));
		JPanel buttonPanel = new JPanel();
		buttonPanel.add(closeButton);
		frame.getContentPane().add(buttonPanel, BorderLayout.PAGE_END);

		// Redirect streams to output window:
		OutputStream out = new OutputStream() {
			@Override
			public void write(int b) throws IOException {
				SwingUtilities.invokeLater(() -> {
					if ((char) b == '\r') {
						int from = textArea.getText().lastIndexOf("\n") + 1;
						int to = textArea.getText().length();
						textArea.replaceRange(null, from, to);
					} else
						textArea.append(String.valueOf((char) b));
				});
			}
		};

		System.setOut(new PrintStream(out, true));
		System.setErr(new PrintStream(out, true));

		frame.pack();
		frame.setVisible(true);
	}

	public static String helpMessage = "Indirect transmission analyser - outputs transmissions that are compatible with the tree topology.\n"
			+ "\n"
			+ "Usage: appstore IndirectTransmissionAnalyser [-help] |  logFile [outputFile]\n"
			+ "\n"
			+ "Option                   Description\n"
			+ "--------------------------------------------------------------\n"
			+ "-help                    Display usage info.\n"
			+ "\n"
			+ "If no output file is specified, output is written to a file\n"
			+ "named 'summary.tree'.";

	/**
	 * Print usage info and exit.
	 */
	public static void printUsageAndExit() {
		System.out.println(helpMessage);
		System.exit(0);
	}

	/**
	 * Display error, print usage and exit with error.
	 */
	public static void printUsageAndError(String errMsg) {
		System.err.println(errMsg);
		System.err.println(helpMessage);
		System.exit(1);
	}

	/**
	 * Retrieve ACGAnnotator options from command line.
	 *
	 * @param args    command line arguments
	 * @param options object to populate with options
	 */
	public static void getCLIOptions(String[] args, TntAnalyserOptions options) {
		int i = 0;
		while (args[i].startsWith("-")) {
			switch (args[i]) {
			case "-help":
				printUsageAndExit();
				break;
			
			default:
				printUsageAndError("Unrecognised command line option '" + args[i] + "'.");
			}
			i += 1;
		}

		if (i >= args.length)
			printUsageAndError("No input file specified.");
		else
			options.inTransmissionTreeFile = new File(args[i]);

		if (i + 1 < args.length)
			options.outFile = new File(args[i + 1]);
		else
			options.outFile = new File("tntAnalysis.log");
	}

	/**
	 * Main method for TnTOrientator. Sets up GUI if needed then uses the
	 * TnTOrientator constructor to actually perform the analysis.
	 *
	 * @param args command line arguments
	 */
	public static void main(String[] args) {
		TntAnalyserOptions options = new TntAnalyserOptions();

		if (args.length == 0) {
			// Retrieve options from GUI:

			try {
				UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
			} catch (ClassNotFoundException | InstantiationException | UnsupportedLookAndFeelException
					| IllegalAccessException e) {
				Log.warning.println("Error setting cross-platform look and feel.");
			}


		} else {
			getCLIOptions(args, options);
		}

		// Run ACGAnnotator
		try {
			new IndirectTransmissionAnalyser2(options);
		} catch (Exception e) {
			if (args.length == 0) {
				JOptionPane.showMessageDialog(null, e.getMessage(),
						"Error", JOptionPane.ERROR_MESSAGE);
			} else {
				System.err.println("Error: " + e.getMessage());
				e.printStackTrace();
				System.err.println();
				System.err.println(helpMessage);
			}

			System.exit(1);
		}
	}

}
