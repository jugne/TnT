package tnt.tntAnnotator;

import java.awt.BorderLayout;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.WindowConstants;

import beast.app.treeannotator.TreeAnnotator;
import beast.core.util.Log;
import beast.evolution.tree.CladeSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeUtils;

public class TransmissionAnalyser extends TreeAnnotator {

	private static class TntAnalyserOptions extends TreeAnnotator {
		File inTransmissionTreeFile;
		List<String> leafIds;
		Double minHeight;
		Double maxHeight;
		File outFile;

		@Override
		public String toString() {
			return "Active options:\n" +
					"Input transmission tree file: " + inTransmissionTreeFile + "\n" +
					"Output file: " + outFile + "\n" +
					"Input sample IDs: " + leafIds + "\n" +
					"Minimum transmission height: " + minHeight + "\n" +
					"Maximum transmission height: " + maxHeight + "\n";

		}

		private String outputFileName() {
			String name = "";
			for (String s : leafIds) {
				name = name + s + "_";
			}

			return name + "tntAnalysis.log";
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

	public TransmissionAnalyser(TntAnalyserOptions options) throws IOException {
 	
    	
        // Display options:
        System.out.println(options + "\n");

        // Initialise reader

		TreeSet trTreeSet = new MemoryFriendlyTreeSet(options.inTransmissionTreeFile.toString(), 0);// new
																									// FastTreeSet(options.inTransmissionTreeFile.toString(),
																									// 0);
		Set<String> leafNodes = new HashSet<String>(options.leafIds);

		PrintStream ps = new PrintStream(options.outFile);
		ps.print("mrca_height");
		for (int i=0; i<options.leafIds.size(); i++) {
			ps.print("_"+options.leafIds.get(i));
		}
		ps.print("\t");

		ps.print("within_interval" + "\t");

		ps.print("before_lower_bound" + "\t");

		for (int i=0; i<options.leafIds.size(); i++) {
			ps.print(options.leafIds.get(i) + "_recipient_mrca" + "\t");
		}

		for (int i = 0; i < options.leafIds.size(); i++) {
			ps.print(options.leafIds.get(i) + "_recipient" + "\t");
		}
		
		for (int i = 0; i < options.leafIds.size(); i++) {
			ps.print(options.leafIds.get(i) + "_n_intermediate_transmissions" + "\t");
		}

		for (int i = 0; i < options.leafIds.size(); i++) {
			ps.print(options.leafIds.get(i) + "_taxa_intermediate_transmissions" + "\t");
		}

		ps.print("other_nodes_in_clade");
		
		ps.print("\n");
		

		trTreeSet.reset();
		while (trTreeSet.hasNext()) {
			Tree tree = trTreeSet.next();
			Node mrca = TreeUtils.getCommonAncestorNode(tree, leafNodes);
			CladeSet clades = new CladeSet(tree);
			int idx = Integer.MIN_VALUE;
//			int i = 0;
			for (int i = 0; i < clades.getCladeCount(); i++) {
				if (mrca.getHeight() == clades.getMeanNodeHeight(i)) {
					idx = i;
					break;
				}
			}
			BitSet clade = clades.get(idx);
			Boolean above_lower = mrca.getHeight() >= options.minHeight;
			Boolean inInterval = above_lower && mrca.getHeight() <= options.maxHeight;
			int inIntervalInt = inInterval ? 1 : 0;
			int aboveLowerBound = above_lower ? 1 : 0;
			ps.print(mrca.getHeight() + "\t" +
					inIntervalInt + "\t" +
					aboveLowerBound + "\t");


			Integer[] isRecipient_mrca = new Integer[options.leafIds.size()];
			Node leftChild = mrca.getChild(0).metaDataString.contains("orientation=donor") ? mrca.getChild(0)
					: mrca.getChild(1);
			for (int j = 0; j < options.leafIds.size(); j++) {
				Set<String> leftDescendants = TreeUtils.getDescendantLeaves(tree, leftChild);
				isRecipient_mrca[j] = leftDescendants.contains(options.leafIds.get(j)) ? 0 : 1;
			}

			for (int i = 0; i < isRecipient_mrca.length; i++) {
				ps.print(isRecipient_mrca[i] + "\t");
			}


			String other_nodes = "";
			Integer[] isRecipient = new Integer[options.leafIds.size()];
			HashMap<Integer, List<String>> separatingBifurcations = new HashMap<>();
			HashMap<Integer, Integer> separatingBifurcationsSize = new HashMap<>();
			List<String> otherNodes_ = new ArrayList<String>();
			for (int n = 0; n < tree.getLeafNodeCount(); n++) {
				for (int j=0; j< options.leafIds.size(); j++) {
					if (tree.getNode(n).getID().equalsIgnoreCase(options.leafIds.get(j))) {
						SeparatingBifurcations s = GetSeparatingBifurcations(tree, mrca, tree.getNode(n));
						separatingBifurcations.put(j, s.taxaList);
						separatingBifurcationsSize.put(j, s.n_transmissions);
						isRecipient[j] = tree.getNode(n).metaDataString.contains("orientation=recipient") ? 1 : 0;
					}

				}
				if (clade.get(n) && !options.leafIds.contains(tree.getNode(n).getID())) {
					otherNodes_.add(tree.getNode(n).getID());
				}
			}

			for (int i = 0; i < isRecipient.length; i++) {
				ps.print(isRecipient[i] + "\t");
			}
			
			for (int i = 0; i < options.leafIds.size(); i++) {
				ps.print(separatingBifurcationsSize.get(i) + "\t");
			}
			
			for (int i = 0; i < options.leafIds.size(); i++) {
				ps.print(separatingBifurcations.get(i) + "\t");
			}

			other_nodes = String.join(",", otherNodes_);
			ps.print(other_nodes);
			ps.print("\n");


        }
		System.out.println("\nDone!");
 

        }

	private SeparatingBifurcations GetSeparatingBifurcations(Tree tree, Node subtreeRoot, Node leaf) {

		Node parent = leaf.getParent();
		Node child = leaf;
		List<String> ans = new ArrayList<>();
		int count = 0;
		while (subtreeRoot.getNr() != parent.getNr()) {
			if (!parent.isFake()) {
				count += 1;
				Node otherChild = parent.getChild(0).getNr() == child.getNr() ? parent.getChild(1) : parent.getChild(0);
				ans.addAll(TreeUtils.getDescendantLeavesSortedList(tree, otherChild));
//				if (otherChild.isFake()) {
//					 Node taxa =  otherChild.getChild(0).getHeight() == otherChild.getHeight() ? otherChild.getChild(0) : otherChild.getChild(1);
//					 ans.add(taxa.getID());
//					 
//				} else if(otherChild.isLeaf()) {
//					ans.add(otherChild.getID());
//				} else {
//					while ()
//				}
					
			}

			
			parent = parent.getParent();
		}
		

		return new SeparatingBifurcations(ans, count);
	}



//	/**
//	 * Use a GUI to retrieve ACGAnnotator options.
//	 *
//	 * @param options options object to populate using GUI
//	 * @return true if options successfully collected, false otherwise
//	 */
//	private static boolean getOptionsGUI(TntAnalyserOptions options) {
//
//		boolean[] canceled = { false };
//
//		JDialog dialog = new JDialog((JDialog) null, true);
//		dialog.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
//		dialog.setLocationRelativeTo(null);
//		dialog.setTitle("Isolation with Migration Annotator");
//
//		JLabel logFileLabel = new JLabel("Isolation with migration species tree file:");
//		JLabel outFileLabel = new JLabel("Output file:");
//
//		JTextField inFilename = new JTextField(20);
//		inFilename.setEditable(false);
//		JButton inFileButton = new JButton("Choose File");
//
//		JTextField outFilename = new JTextField(20);
//		outFilename.setText(options.outFile.getName());
//		outFilename.setEditable(false);
//		JButton outFileButton = new JButton("Choose File");
//
//		Container cp = dialog.getContentPane();
//		BoxLayout boxLayout = new BoxLayout(cp, BoxLayout.PAGE_AXIS);
//		cp.setLayout(boxLayout);
//
//		JPanel mainPanel = new JPanel();
//
//		GroupLayout layout = new GroupLayout(mainPanel);
//		mainPanel.setLayout(layout);
//		layout.setAutoCreateGaps(true);
//		layout.setAutoCreateContainerGaps(true);
//
//		layout.setHorizontalGroup(layout.createSequentialGroup()
//				.addGroup(layout.createParallelGroup()
//						.addComponent(logFileLabel)
//						.addComponent(outFileLabel))
//				.addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
//						.addComponent(inFilename)
//						.addComponent(outFilename))
//				.addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
//						.addComponent(inFileButton)
//						.addComponent(outFileButton)));
//
//		layout.setVerticalGroup(layout.createSequentialGroup()
//				.addGroup(layout.createParallelGroup()
//						.addComponent(logFileLabel)
//						.addComponent(inFilename,
//								GroupLayout.PREFERRED_SIZE,
//								GroupLayout.DEFAULT_SIZE,
//								GroupLayout.PREFERRED_SIZE)
//						.addComponent(inFileButton))
//				.addGroup(layout.createParallelGroup()
//						.addComponent(outFileLabel)
//						.addComponent(outFilename,
//								GroupLayout.PREFERRED_SIZE,
//								GroupLayout.DEFAULT_SIZE,
//								GroupLayout.PREFERRED_SIZE)
//						.addComponent(outFileButton)));
//
//		mainPanel.setBorder(new EtchedBorder());
//		cp.add(mainPanel);
//
//		JPanel buttonPanel = new JPanel();
//
//		JButton runButton = new JButton("Analyze");
//		runButton.addActionListener((e) -> {
//			dialog.setVisible(false);
//		});
//		runButton.setEnabled(false);
//		buttonPanel.add(runButton);
//
//		JButton cancelButton = new JButton("Quit");
//		cancelButton.addActionListener((e) -> {
//			dialog.setVisible(false);
//			canceled[0] = true;
//		});
//		buttonPanel.add(cancelButton);
//
//		JFileChooser inFileChooser = new JFileChooser();
//		inFileButton.addActionListener(e -> {
//			inFileChooser.setDialogTitle("Select Transmission Trees to orientate");
//			if (options.inFile == null)
//				inFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
//			int returnVal = inFileChooser.showOpenDialog(dialog);
//
//			if (returnVal == JFileChooser.APPROVE_OPTION) {
//				options.inFile = inFileChooser.getSelectedFile();
//				inFilename.setText(inFileChooser.getSelectedFile().getName());
//				runButton.setEnabled(true);
//			}
//		});
//
//		JFileChooser outFileChooser = new JFileChooser();
//		outFileButton.addActionListener(e -> {
//			outFileChooser.setDialogTitle("Select output file name.");
//			if (options.inFile != null)
//				outFileChooser.setCurrentDirectory(options.inFile);
//			else
//				outFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
//
//			outFileChooser.setSelectedFile(options.outFile);
//			int returnVal = outFileChooser.showOpenDialog(dialog);
//
//			if (returnVal == JFileChooser.APPROVE_OPTION) {
//				options.outFile = outFileChooser.getSelectedFile();
//				outFilename.setText(outFileChooser.getSelectedFile().getName());
//			}
//		});
//
//		cp.add(buttonPanel);
//
//		dialog.pack();
//		dialog.setResizable(false);
//		dialog.setVisible(true);
//
//		return !canceled[0];
//	}

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

	public static String helpMessage = "Transmission Tree Orientator - orientates transmission from TNT package tree based on the metadata.\n"
			+ "\n"
			+ "Usage: appstore ACGAnnotator [-help | [options] logFile [outputFile]\n"
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
			case "-leafIds":
				options.leafIds = new ArrayList<>(Arrays.asList(args[i + 1].split(",")));
				i += 1;
				break;
			case "-minHeight":
				options.minHeight = new Double(args[i + 1]);
				i += 1;
				break;	
			case "-maxHeight":
				options.maxHeight = new Double(args[i + 1]);
				i += 1;
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
			options.outFile = new File(options.outputFileName());
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

//			try {
//				SwingUtilities.invokeAndWait(() -> {
//					if (!getOptionsGUI(options))
//						System.exit(0);
//
//					setupGUIOutput();
//				});
//			} catch (InterruptedException | InvocationTargetException e) {
//				e.printStackTrace();
//			}

		} else {
			getCLIOptions(args, options);
		}

		// Run ACGAnnotator
		try {
			new TransmissionAnalyser(options);
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
