package tnt.orientation;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.evolution.tree.Node;
import distribution.GeneTree;
import starbeast2.SpeciesTreeInterface;

@Description("Based on the SpeciesTreeLogger class, but without node sorting")
public class OrientationTest extends BEASTObject implements Loggable {
	final public Input<SpeciesTreeInterface> speciesTreeInput = new Input<>("transmissionTree",
			"The species tree to be logged.", Validate.REQUIRED);

	final public Input<List<GeneTree>> geneTreeInput = new Input<>("geneTree", "Gene tree within the species tree.",
			new ArrayList<>());

    private DecimalFormat df;

	private HashMap<Integer, Pattern> rx;
	private int[] freq;

    @Override
    public void initAndValidate() {
    }

    @Override
    public void init(PrintStream out) {
        SpeciesTreeInterface speciesTree = speciesTreeInput.get();
        speciesTree.init(out);

		freq = new int[18];

		String number = "\\:(\\d+(\\.\\d+)?)(E\\-\\d+)?";

		rx = new HashMap<Integer, Pattern>();

		// Group 1. Total probability for non-oriented topology ((C,B),A): 77.8327
		// (A,(B,C))
		rx.put(1, Pattern
				.compile("\\(A(?!(\\:0.0\\,))" + number + "\\,\\(B(?!(\\:0.0\\,))" + number + "\\,C" + number
						+ "\\)(.*)\\)(.*)"));

		// (A,(C,B))
		rx.put(2, Pattern
				.compile("\\(A(?!(\\:0.0\\,))" + number + "\\,\\(C" + number + "\\,B(?!(\\:0.0\\)))" + number
						+ "\\)(.*)\\)(.*)"));

		// ((B,C),A)
		rx.put(3,
				Pattern.compile(
						"\\(\\(B" + number + "\\,C" + number + "\\)(.*)\\,A(?!(\\:0.0\\)))" + number + "\\)(.*)"));

		// ((C,B),A)
		rx.put(4, Pattern
				.compile("\\(\\(C" + number + "\\,B" + number + "\\)(.*)\\,A(?!(\\:0.0\\)))" + number + "\\)(.*)"));


		// Group 2. Total probability for non-oriented topology ((C,A),B): 4.3189
		// (B, (A,C))
		rx.put(5, Pattern.compile("\\(B" + number + "\\,\\(A" + number + "\\,C" + number + "\\)(.*)\\)(.*)"));
		// (B,(C,A))
		rx.put(6, Pattern.compile("\\(B" + number + "\\,\\(C" + number + "\\,A" + number + "\\)(.*)\\)(.*)"));
		// ((A,C),B)
		rx.put(7, Pattern.compile("\\(\\(A" + number + "\\,C" + number + "\\)(.*)\\,B" + number + "\\)(.*)"));
		// ((C,A),B)
		rx.put(8, Pattern.compile("\\(\\(C" + number + "\\,A" + number + "\\)(.*)\\,B" + number + "\\)(.*)"));

		// Group 3. Total probability for non-oriented topology (C,(B,A)): 4.3189
		// (C,(A,B))
		rx.put(9, Pattern
				.compile("\\(C" + number + "\\,\\(A(?!\\:0.0\\,))" + number + "\\,B" + number + "\\)(.*)\\)(.*)"));
		// (C,(B,A))
		rx.put(10, Pattern
				.compile("\\(C" + number + "\\,\\(B" + number + "\\,\\(A(?!\\:0.0\\,))" + number + "\\)(.*)\\)(.*)"));
		// ((B,A), C)
		rx.put(11, Pattern
				.compile("\\(\\(B" + number + "\\,\\(A(?!\\:0.0\\,))" + number + "\\)(.*)\\,C" + number + "\\)(.*)"));
		// ((A,B), C)
		rx.put(12, Pattern
				.compile("\\(\\(A(?!\\:0.0\\,))" + number + "\\,B" + number + "\\)(.*)\\,C" + number + "\\)(.*)"));

		// Group 4. Total probability for non-oriented topology ((B,C)A): 7.8642
		// ((B,C)A)
		rx.put(13, Pattern.compile("\\(\\(B" + number + "\\,C" + number + "\\)(.*)\\,A:0.0\\)(.*)"));
		// ((C,B)A)
		rx.put(14, Pattern.compile("\\(\\(C" + number + "\\,B" + number + "\\)(.*)\\,A:0.0\\)(.*)"));
		// (A(B,C))
		rx.put(15, Pattern.compile("\\(A:0.0\\,\\(B" + number + "\\,C" + number + "\\)(.*)\\)(.*)"));
		// (A(C,B))
		rx.put(16, Pattern.compile("\\(A:0.0\\,\\(C" + number + "\\,B" + number + "\\)(.*)\\)(.*)"));

		// Group 5. Total probability for non-oriented topology ((C)B,A): 0.6930
		// (A,B(C))
		rx.put(17, Pattern.compile("\\(A(?!(\\:0.0\\,))" + number + "\\,\\(B:0.0\\,C" + number + "\\)(.*)\\)(.*)"));
		// (A,(C)B)
		rx.put(18, Pattern.compile("\\(A(?!(\\:0.0\\,))" + number + "\\,\\(C" + number + "\\,B:0.0\\)(.*)\\)(.*)"));
		// (B(C),A)
		rx.put(19, Pattern
				.compile("\\(\\(B:0.0\\,C" + number + "\\)(.*)\\,,A(?!(\\\\:0.0\\\\,))\" + number +\\)(.*)"));
		// ((C)B,A)
		rx.put(20, Pattern.compile("\\(\\(C" + number + "\\,B:0.0\\)(.*)\\,A(?!(\\:0.0\\,))" + number + "\\)(.*)"));


		// Group 6. Total probability of non-oriented topology (C,(B)A): 0.6930
		// (C,A(B))
		rx.put(21, Pattern.compile("\\(C(?!(\\:0.0\\,))" + number + "\\,\\(A:0.0\\,B" + number + "\\)(.*)\\)(.*)"));
		// (C,(B)A)
		rx.put(22, Pattern.compile("\\(C(?!(\\:0.0\\,))" + number + "\\,\\(B" + number + "\\,A:0.0\\)(.*)\\)(.*)"));
		// (A(B), C)
		rx.put(23, Pattern.compile("\\(\\(A:0.0\\,B" + number + "\\)\\, C(?!(\\:0.0\\,))" + number + "\\)(.*)\\)(.*)"));
		// ((B)A,C)
		rx.put(24, Pattern.compile("\\(C(?!(\\\\:0.0\\\\,))" + number + "\\,\\(A:0.0\\,B" + number + "\\)(.*)\\)(.*)"));

		// Group 7. Total probability of non-oriented topology (((C)B)A): 0.4135

    }

    @Override
    public void log(long nSample, PrintStream out) {
        // make sure we get the current version of the inputs
        SpeciesTreeInterface speciesTree = speciesTreeInput.get();
        SpeciesTreeInterface tree = (SpeciesTreeInterface) speciesTree.getCurrent();

        // write out the log tree with meta data
        out.print("tree STATE_" + nSample + " = ");
//        tree.getRoot().sort();
		out.print(toNewick(tree.getRoot()));
        //out.print(tree.getRoot().toShortNewick(false));
        out.print(";");

		String newick = toNewick(tree.getRoot());

		int[] duplicate = Arrays.copyOf(freq, freq.length);
		for (int i = 1; i <= freq.length; i++) {
			Matcher m = rx.get(i).matcher(newick);
			if (m.matches()) {
				freq[i - 1] += 1;
				if (i == 17) {
					System.out.println(newick);
				}
				if (i == 18) {
					System.out.println(newick);
				}
			}

		}

		if (Arrays.equals(duplicate, freq)) {
			System.out.println("Not Assigned: " + newick);
			System.exit(0);
		}

		if (Arrays.stream(freq).sum() != Arrays.stream(duplicate).sum() + 1) {
			System.out.println("Tree assigned to more than one topology");
			System.exit(0);
		}

		System.out.println(Arrays.toString(freq));

    }

    /**
     * Appends a double to the given StringBuffer, formatting it using
     * the private DecimalFormat instance, if the input 'dp' has been
     * given a non-negative integer, otherwise just uses default
     * formatting.
     * @param buf
     * @param d
     */
    private void appendDouble(StringBuffer buf, double d) {
        if (df == null) {
            buf.append(d);
        } else {
            buf.append(df.format(d));
        }
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
        } else {
//            buf.append(node.getNr() + 1);
        }

		if (node.getID() == null) {
			buf.append(node.getNr() + 1);
		} else {
			buf.append(node.getID());
		}
        buf.append(":");

        double nodeLength;
        if (node.isRoot()) {
            nodeLength = getTreeHeight() - node.getHeight();
        } else {
            nodeLength = node.getLength();
        }
		appendDouble(buf, nodeLength);

        return buf.toString();
    }

    // uses the height of the tallest species or gene tree
    private double getTreeHeight() {
        double speciesTreeHeight = speciesTreeInput.get().getRoot().getHeight();

        for (GeneTree gt: geneTreeInput.get()) {
            speciesTreeHeight = Double.max(speciesTreeHeight, gt.getRoot().getHeight());
        }

        return speciesTreeHeight;
    }

    @Override
    public void close(PrintStream out) {
        SpeciesTreeInterface speciesTree = speciesTreeInput.get();
        speciesTree.close(out);
    }

}

