package tnt.mapping;

import bdmmprime.parameterization.SkylineVectorParameter;
import bdmmprime.parameterization.TimedParameter;
import bdmmprime.parameterization.TypeSet;
import beast.app.BeastMCMC;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import feast.fileio.logfileiterator.LogFileIterator;
import feast.fileio.logfileiterator.LogFileRealParameter;
import feast.fileio.logfileiterator.TraceLogFileState;
import feast.fileio.logfileiterator.TreeLogFileState;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import tnt.transmissionTree.TransmissionTree;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import java.awt.*;
import java.io.*;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

//TODO make work with different parameterization (now only Epi)
//TODO make work with non-present rho sampling
/**
 * @author ugne.stolz@protonmail.com
 * @date 11.04.22
 */
public class ReMapTool {


    private static class ReMapOptions {
        File xmlFile;
        File logFile;
        File treeFile;
        File outFile;

        @Override
        public String toString() {
            return "Active options:\n" +
                    "TnT XML file: " + xmlFile + "\n" +
                    "Trees file: " + treeFile + "\n" +
                    "Log file: " + logFile + "\n" +
                    "Output file: " + outFile;
        }
    }

    public ReMapTool(ReMapOptions options) throws Exception{

        // Display options:
        System.out.println(options + "\n");

        BufferedReader reader = new BufferedReader(new FileReader(options.logFile));
        String header = reader.readLine();
        reader.close();

        String[] headerSplit = header.split("\\t");
        Pattern originPattern = Pattern.compile("origin");
        Pattern R0Pattern = Pattern.compile("R0");
        Pattern deltaPattern = Pattern.compile("becominUninfectiousRate");
        Pattern sPattern = Pattern.compile("samplingProportion");
        Pattern rPattern = Pattern.compile("removalProb");
        Pattern rhoPattern = Pattern.compile("rhoSampling");

        Boolean rhoUsed = true;


        // get header strings for all logger types
        String origin = null;
        String R0 =null;
        String delta = null;
        String s =null;
        String r =null;
        String rho = null;

        for (String param : headerSplit) {
            Matcher originMatcher = originPattern.matcher(param);
            Matcher R0Matcher = R0Pattern.matcher(param);
            Matcher deltaMatcher = deltaPattern.matcher(param);
            Matcher sMatcher = sPattern.matcher(param);
            Matcher rMatcher = rPattern.matcher(param);
            Matcher rhoMatcher = rhoPattern.matcher(param);

            if (originMatcher.find())
                origin = originMatcher.group(0);
            if (R0Matcher.find())
                R0 = R0Matcher.group(0);
            if (deltaMatcher.find())
                delta = deltaMatcher.group(0);
            if (sMatcher.find())
                s = sMatcher.group(0);
            if (rMatcher.find())
                r = rMatcher.group(0);
            if (rhoMatcher.find())
                rho = rhoMatcher.group(0);
        }

        if (origin == null)
            throw new Exception("Origin not found in log file.");
        if (R0 == null)
            throw new Exception("R0 not found in log file.");
        if (delta == null)
            throw new Exception("Become uninfectious rate not found in log file.");
        if (s == null)
            throw new Exception("Sampling proportion not found in log file.");
        if (r == null)
            throw new Exception("Removal probability not found in log file.");
        if (rho == null)
            rhoUsed = false;


        DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
        DocumentBuilder docBuilder = docFactory.newDocumentBuilder();
        Document tntXML = docBuilder.parse(options.xmlFile);

        NodeList taxonSets = tntXML.getElementsByTagName("taxonsets");
        List<Element> taxonsets = new ArrayList<>();
        for(int l=0; l < taxonSets.getLength(); l++){
            taxonsets.add((Element) taxonSets.item(l));
        }

        Element offset = (Element) tntXML.getElementsByTagName("finalSampleOffset").item(0);

        // retrieve the element 'run'
        Element run = (Element) tntXML.getElementsByTagName("run").item(0);

        // remove some elements
        Element birthDeathSkylineModel = (Element) tntXML.getElementsByTagName("BirthDeathSkylineModel").item(0);
        Element parameterization = (Element) tntXML.getElementsByTagName("parameterization").item(0);

        // get parent node of run which will be used later to create a new run node
        Element beast = (Element) run.getParentNode();

        // remove the old run node
        beast.removeChild(run);
        beast.removeChild(birthDeathSkylineModel);
        beast.removeChild(parameterization);

        Element feastRun = tntXML.createElement("run");
        beast.appendChild(feastRun);
        feastRun.setAttribute("spec", LogFileIterator.class.getCanonicalName());

        // trace log state
        Element traceLogState = tntXML.createElement("logFileState");
        traceLogState.setAttribute("spec", TraceLogFileState.class.getCanonicalName());
        traceLogState.setAttribute("logFileName", options.logFile.getName());
        feastRun.appendChild(traceLogState);

        Element originLogFileEntries = tntXML.createElement("logFileEntry");
        originLogFileEntries.setAttribute("spec", LogFileRealParameter.class.getCanonicalName());
        originLogFileEntries.setAttribute("fieldName", origin);
        Element originFieldParameters = tntXML.createElement("fieldParameter");
        originFieldParameters.setAttribute("id", origin);
        originFieldParameters.setAttribute("spec", RealParameter.class.getCanonicalName());
        originFieldParameters.setAttribute("value", "0.0");

        originLogFileEntries.appendChild(originFieldParameters);
        traceLogState.appendChild(originLogFileEntries);

        Element R0LogFileEntries = tntXML.createElement("logFileEntry");
        R0LogFileEntries.setAttribute("spec", LogFileRealParameter.class.getCanonicalName());
        R0LogFileEntries.setAttribute("fieldName", R0);
        Element R0FieldParameters = tntXML.createElement("fieldParameter");
        R0FieldParameters.setAttribute("id", R0);
        R0FieldParameters.setAttribute("spec", RealParameter.class.getCanonicalName());
        R0FieldParameters.setAttribute("value", "0.0");

        R0LogFileEntries.appendChild(R0FieldParameters);
        traceLogState.appendChild(R0LogFileEntries);

        Element deltaLogFileEntries = tntXML.createElement("logFileEntry");
        deltaLogFileEntries.setAttribute("spec", LogFileRealParameter.class.getCanonicalName());
        deltaLogFileEntries.setAttribute("fieldName", delta);
        Element deltaFieldParameters = tntXML.createElement("fieldParameter");
        deltaFieldParameters.setAttribute("id", delta);
        deltaFieldParameters.setAttribute("spec", RealParameter.class.getCanonicalName());
        deltaFieldParameters.setAttribute("value", "0.0");

        deltaLogFileEntries.appendChild(deltaFieldParameters);
        traceLogState.appendChild(deltaLogFileEntries);

        Element sLogFileEntries = tntXML.createElement("logFileEntry");
        sLogFileEntries.setAttribute("spec", LogFileRealParameter.class.getCanonicalName());
        sLogFileEntries.setAttribute("fieldName", s);
        Element sFieldParameters = tntXML.createElement("fieldParameter");
        sFieldParameters.setAttribute("id", s);
        sFieldParameters.setAttribute("spec", RealParameter.class.getCanonicalName());
        sFieldParameters.setAttribute("value", "0.0");

        sLogFileEntries.appendChild(sFieldParameters);
        traceLogState.appendChild(sLogFileEntries);

        Element rLogFileEntries = tntXML.createElement("logFileEntry");
        rLogFileEntries.setAttribute("spec", LogFileRealParameter.class.getCanonicalName());
        rLogFileEntries.setAttribute("fieldName", r);
        Element rFieldParameters = tntXML.createElement("fieldParameter");
        rFieldParameters.setAttribute("id", r);
        rFieldParameters.setAttribute("spec", RealParameter.class.getCanonicalName());
        rFieldParameters.setAttribute("value", "0.0");

        rLogFileEntries.appendChild(rFieldParameters);
        traceLogState.appendChild(rLogFileEntries);

        if (rhoUsed){
            Element rhoLogFileEntries = tntXML.createElement("logFileEntry");
            rhoLogFileEntries.setAttribute("spec", LogFileRealParameter.class.getCanonicalName());
            rhoLogFileEntries.setAttribute("fieldName", rho);
            Element rhoFieldParameters = tntXML.createElement("fieldParameter");
            rhoFieldParameters.setAttribute("id", rho);
            rhoFieldParameters.setAttribute("spec", RealParameter.class.getCanonicalName());
            rhoFieldParameters.setAttribute("value", "0.0");

            rhoLogFileEntries.appendChild(rhoFieldParameters);
            traceLogState.appendChild(rhoLogFileEntries);
        }
        // network log state
        Element treeLogState = tntXML.createElement("logFileState");
        treeLogState.setAttribute("spec", TreeLogFileState.class.getCanonicalName());
        treeLogState.setAttribute("logFileName", options.treeFile.getName());

        Element tree = tntXML.createElement("tree");
        tree.setAttribute("spec", TransmissionTree.class.getCanonicalName());
        tree.setAttribute("id", "tree");
        tree.setAttribute("taxonset", "@taxonsuperset");
        tree.setAttribute("adjustTreeNodeHeights", "false");
        tree.setAttribute("trait", "@tipDates");
        treeLogState.appendChild(tree);

//        Element taxonSet = tntXML.createElement("readTaxonSet");
//        taxonSet.setAttribute("spec", Boolean.class.getCanonicalName());
//        taxonSet.setAttribute("value", "0.0");
//        treeLogState.appendChild(taxonSet);

        feastRun.appendChild(treeLogState);

        // loggers

        // network logger
        Element treeLogger = tntXML.createElement("logger");
        treeLogger.setAttribute("spec", "Logger");
        treeLogger.setAttribute("logEvery", "1");
        treeLogger.setAttribute("mode", "tree");
        treeLogger.setAttribute("fileName", "$(filebase).trees");

        Element treeLog = tntXML.createElement("log");
        treeLog.setAttribute("spec", MappedTree.class.getCanonicalName());
        treeLog.setAttribute("tree", "@tree");
        treeLog.setAttribute("mapOnInit", "false");
        treeLog.setAttribute("id", "mappedTree");
        for (Element e : taxonsets){
            treeLog.appendChild(e);
        }

        Element param = tntXML.createElement("parameterization");
        param.setAttribute("id", "parameterization");
        param.setAttribute("spec", "EpiParameterization");

//        <typeSet id="typeSet" spec="bdmmprime.parameterization.TypeSet" value="0"/>
        Element typeSet = tntXML.createElement("typeSet");
        typeSet.setAttribute("id","typeSet");
        typeSet.setAttribute("spec", TypeSet.class.getCanonicalName());
        typeSet.setAttribute("value", "0");

        param.appendChild(typeSet);

//        <origin id="origin" spec="RealParameter" value="30." lower="0."/>
        Element originLog = tntXML.createElement("origin");
        originLog.setAttribute("idref",origin);
        param.appendChild(originLog);

//        <R0 spec="SkylineVectorParameter" typeSet="@typeSet">
//      		<skylineValues id="R0" spec="RealParameter" value="3.33180010714568"/>
//      	</R0>
        Element R0Log = tntXML.createElement("R0");
        R0Log.setAttribute("spec", SkylineVectorParameter.class.getCanonicalName());
        R0Log.setAttribute("typeSet", "@typeSet");
        Element R0LogVal = tntXML.createElement("skylineValues");
        R0LogVal.setAttribute("idref", R0);
        R0Log.appendChild(R0LogVal);
        param.appendChild(R0Log);

//        <becomeUninfectiousRate spec="SkylineVectorParameter" typeSet="@typeSet">
//      		<skylineValues id="becominUninfectiousRate" spec="RealParameter" value="1.09589617969468"/>
//      	</becomeUninfectiousRate>
        Element deltaLog = tntXML.createElement("becomeUninfectiousRate");
        deltaLog.setAttribute("spec", SkylineVectorParameter.class.getCanonicalName());
        deltaLog.setAttribute("typeSet", "@typeSet");
        Element deltaLogVal = tntXML.createElement("skylineValues");
        deltaLogVal.setAttribute("idref", delta);
        deltaLog.appendChild(deltaLogVal);
        param.appendChild(deltaLog);

//        <samplingProportion spec="SkylineVectorParameter" typeSet="@typeSet">
//      		<skylineValues id="samplingProportion" spec="RealParameter" value="0.752627511392348"/>
//      	</samplingProportion>
        Element sLog = tntXML.createElement("samplingProportion");
        sLog.setAttribute("spec", SkylineVectorParameter.class.getCanonicalName());
        sLog.setAttribute("typeSet", "@typeSet");
        Element sLogVal = tntXML.createElement("skylineValues");
        sLogVal.setAttribute("idref", s);
        sLog.appendChild(sLogVal);
        param.appendChild(sLog);

//        <removalProb spec="SkylineVectorParameter" typeSet="@typeSet">
//      		<skylineValues id="removalProb" spec="RealParameter" value="0.226106221601367"/>
//      	</removalProb>
        Element rLog = tntXML.createElement("removalProb");
        rLog.setAttribute("spec", SkylineVectorParameter.class.getCanonicalName());
        rLog.setAttribute("typeSet", "@typeSet");
        Element rLogVal = tntXML.createElement("skylineValues");
        rLogVal.setAttribute("idref", r);
        rLog.appendChild(rLogVal);
        param.appendChild(rLog);

        if (rhoUsed) {
//        <rhoSampling spec="TimedParameter" typeSet="@typeSet" timesAreAges="true" origin="@origin">
//        	<times spec="RealParameter" value="0.0"/>
//        	<values id="rho" spec="RealParameter" value="0"/>
//      	</rhoSampling>
            Element rhoLog = tntXML.createElement("rhoSampling");
            rhoLog.setAttribute("spec", TimedParameter.class.getCanonicalName());
            rhoLog.setAttribute("typeSet", "@typeSet");
            rhoLog.setAttribute("timesAreAges", "true");
            rhoLog.setAttribute("origin", "@origin");
            Element rhoLogTimes = tntXML.createElement("times");
            rhoLogTimes.setAttribute("spec", RealParameter.class.getCanonicalName());
            rhoLogTimes.setAttribute("value", "0.0");
            Element rhoLogVal = tntXML.createElement("values");
            rhoLogVal.setAttribute("idref", rho);
            rhoLog.appendChild(rhoLogTimes);
            rhoLog.appendChild(rLogVal);
            param.appendChild(rhoLog);
        }

        treeLog.appendChild(param);
        treeLog.appendChild(offset);

        Element nHiddenTransmissions = tntXML.createElement("nHiddenEvents");
        nHiddenTransmissions.setAttribute("id", "nHiddenEvents");
        nHiddenTransmissions.setAttribute("spec", IntegerParameter.class.getCanonicalName());
        nHiddenTransmissions.setAttribute("value", "0");
        treeLog.appendChild(nHiddenTransmissions);

        treeLogger.appendChild(treeLog);
        feastRun.appendChild(treeLogger);

        // stat logger
        Element mapStat = tntXML.createElement("logger");
        mapStat.setAttribute("spec", "Logger");
        mapStat.setAttribute("logEvery", "1");
        mapStat.setAttribute("fileName", "$(filebase).stats.log");

        Element hidStatLog = tntXML.createElement("log");
        hidStatLog.setAttribute("idref", "nHiddenEvents");
        mapStat.appendChild(hidStatLog);
        feastRun.appendChild(mapStat);

        // create the xml file
        // transform the DOM Object to an XML File
        TransformerFactory transformerFactory = TransformerFactory.newInstance();
        Transformer transformer = transformerFactory.newTransformer();
        transformer.setOutputProperty(OutputKeys.INDENT, "yes");
        transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "4");
        DOMSource domSource = new DOMSource(tntXML);
        StreamResult streamResult = new StreamResult(options.outFile);

        transformer.transform(domSource, streamResult);

        System.out.println("Done creating re-map XML File");

        String[] args = new String[2];
        args[0] = "-overwrite";
        args[1] = options.outFile.toString();
        BeastMCMC.main(args);

    }


    /**
     * Use a GUI to retrieve ACGAnnotator options.
     *
     * @param options options object to populate using GUI
     * @return true if options successfully collected, false otherwise
     */
    private static boolean getOptionsGUI(ReMapOptions options) {

        boolean[] canceled = {false};

        JDialog dialog = new JDialog((JDialog)null, true);
        dialog.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        dialog.setLocationRelativeTo(null);
        dialog.setTitle("TnT Re-Mapper");

        JLabel xmlFileLabel = new JLabel("TnT XML file:");
        JLabel netFileLabel = new JLabel("Tree file:");
        JLabel logFileLabel = new JLabel("Log file:");
        JLabel outFileLabel = new JLabel("Output file:");

        JTextField xmlFilename = new JTextField(20);
        xmlFilename.setEditable(false);
        JButton xmlFileButton = new JButton("Choose File");

        JTextField netFilename = new JTextField(20);
        netFilename.setEditable(false);
        JButton netFileButton = new JButton("Choose File");

        JTextField logFilename = new JTextField(20);
        logFilename.setEditable(false);
        JButton logFileButton = new JButton("Choose File");

        JTextField outFilename = new JTextField(20);
        outFilename.setEditable(false);
        JButton outFileButton = new JButton("Choose File");


        Container cp = dialog.getContentPane();
        BoxLayout boxLayout = new BoxLayout(cp, BoxLayout.PAGE_AXIS);
        cp.setLayout(boxLayout);

        JPanel mainPanel = new JPanel();

        GroupLayout layout = new GroupLayout(mainPanel);
        mainPanel.setLayout(layout);
        layout.setAutoCreateGaps(true);
        layout.setAutoCreateContainerGaps(true);

        layout.setHorizontalGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup()
                        .addComponent(xmlFileLabel)
                        .addComponent(netFileLabel)
                        .addComponent(logFileLabel)
                        .addComponent(outFileLabel))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(xmlFilename)
                        .addComponent(netFilename)
                        .addComponent(logFilename)
                        .addComponent(outFilename))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(xmlFileButton)
                        .addComponent(netFileButton)
                        .addComponent(logFileButton)
                        .addComponent(outFileButton))
        );

        layout.setVerticalGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup()
                        .addComponent(xmlFileLabel)
                        .addComponent(xmlFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
                        .addComponent(xmlFileButton))
                .addGroup(layout.createParallelGroup()
                        .addComponent(netFileLabel)
                        .addComponent(netFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
                        .addComponent(netFileButton))
                .addGroup(layout.createParallelGroup()
                        .addComponent(logFileLabel)
                        .addComponent(logFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
                        .addComponent(logFileButton))
                .addGroup(layout.createParallelGroup()
                        .addComponent(outFileLabel)
                        .addComponent(outFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
                        .addComponent(outFileButton))
        );

        mainPanel.setBorder(new EtchedBorder());
        cp.add(mainPanel);

        JPanel buttonPanel = new JPanel();

        JButton runButton = new JButton("Run");
        runButton.addActionListener((e) -> dialog.setVisible(false));
        runButton.setEnabled(false);
        buttonPanel.add(runButton);

        JButton cancelButton = new JButton("Quit");
        cancelButton.addActionListener((e) -> {
            dialog.setVisible(false);
            canceled[0] = true;
        });
        buttonPanel.add(cancelButton);

        JFileChooser xmlFileChooser = new JFileChooser();
        xmlFileButton.addActionListener(e -> {
            xmlFileChooser.setDialogTitle("Select TnT XML file");
            if (options.xmlFile == null)
                xmlFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
            int returnVal = xmlFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                options.xmlFile = xmlFileChooser.getSelectedFile();
                xmlFilename.setText(xmlFileChooser.getSelectedFile().getName());
                runButton.setEnabled(true);
            }
        });

        JFileChooser netFileChooser = new JFileChooser();
        netFileButton.addActionListener(e -> {
            netFileChooser.setDialogTitle("Select TnT transmission tree file to remap");
            if (options.treeFile == null)
                netFileChooser.setCurrentDirectory(options.xmlFile);
            else
                netFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
            int returnVal = netFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                options.treeFile = netFileChooser.getSelectedFile();
                netFilename.setText(netFileChooser.getSelectedFile().getName());
                runButton.setEnabled(true);
            }
        });

        JFileChooser logFileChooser = new JFileChooser();
        logFileButton.addActionListener(e -> {
            logFileChooser.setDialogTitle("Select TnT log file");
            if (options.logFile == null)
                logFileChooser.setCurrentDirectory(options.xmlFile);
            else
                logFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
            int returnVal = logFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                options.logFile = logFileChooser.getSelectedFile();
                logFilename.setText(logFileChooser.getSelectedFile().getName());
                runButton.setEnabled(true);
            }
        });

        JFileChooser outFileChooser = new JFileChooser();
        outFileButton.addActionListener(e -> {
            outFileChooser.setDialogTitle("Select output XML file name");
            outFileChooser.setCurrentDirectory(Objects.requireNonNullElseGet(options.xmlFile, () -> new File(System.getProperty("user.dir"))));

//            outFileChooser.setSelectedFile(options.outFile);
            int returnVal = outFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                options.outFile = outFileChooser.getSelectedFile();
                outFilename.setText(outFileChooser.getSelectedFile().getName());
                runButton.setEnabled(true);
            }
        });

        cp.add(buttonPanel);

        dialog.pack();
        dialog.setResizable(false);
        dialog.setVisible(true);

        return !canceled[0];
    }

    /**
     * Prepare JFrame to which ACGAnnotator output streams will be
     * directed.
     */
    private static void setupGUIOutput() {

        JFrame frame = new JFrame();
        frame.setTitle("TnT Re-Map Tool");
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
            public void write(int b) {
                SwingUtilities.invokeLater(() -> {
                    if ((char)b == '\r') {
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

    public static String helpMessage =
            "ReMapTool - maps hidden transmission events that change the host on a tree branch.\n"
                    + "\n"
                    + "Usage: appstore ReMapTool [-help | [options] logFile [outputFile]\n"
                    + "\n"
                    + "Option                   Description\n"
                    + "--------------------------------------------------------------\n"
                    + "-help                    Display usage info.\n"
                    + "-xml 					XML file of the original TnT analysis.\n"
                    + "-tree				    *.trees transmission trees file produced by TnT.\n"
                    + "-log     				log file produced by TnT.\n"
                    + "-out     				Name of the re-mapping XML file.\n"
                    + "                         Mapped tree's log file naming will follow\n"
                    + "                         output xml naming scheme.\n"
                    + "\n";

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
     * Retrieve ReMapTool options from command line.
     *
     * @param args command line arguments
     * @param options object to populate with options
     */
    public static void getCLIOptions(String[] args, ReMapOptions options) {
        int i=0;
        while (args.length > i && args[i].startsWith("-")) {
            switch(args[i]) {
                case "-help":
                    printUsageAndExit();
                    break;

                case "-xml":
                    if (args.length<=i+1)
                        printUsageAndError("-xml must be followed by a xml file path.");

                    try {
                        options.xmlFile = new File(args[i + 1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("Error parsing XML file.");
                    }

                    i += 1;
                    break;


                case "-tree":
                    if (args.length<=i+1) {
                        printUsageAndError("-tree must be followed by a transmission tree file path.");
                    }

                    try {
                        options.treeFile = new File(args[i + 1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("Error parsing network file.");
                    }

                    i += 1;
                    break;

                case "-log":
                    if (args.length <= i + 1) {
                        printUsageAndError("-log must be followed by a log file path.");
                    }

                    try {
                        options.logFile = new File(args[i + 1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("Error parsing log file.");
                    }

                    i += 1;
                    break;

                case "-out":
                    if (args.length<=i+1) {
                        printUsageAndError("-out must be followed by an output file path.");
                    }

                    try {
                        options.outFile = new File(args[i + 1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("Error parsing output file path.");
                    }

                    i += 1;
                    break;


                default:
                    printUsageAndError("Unrecognised command line option '" + args[i] + "'.");
            }

            i += 1;
        }
    }

    /**
     * Main method for ACGAnnotator.  Sets up GUI if needed then
     * uses the ACGAnnotator constructor to actually perform the analysis.
     *
     * @param args command line arguments
     */
    public static void main(String[] args) {
        ReMapOptions options = new ReMapOptions();

        if (args.length == 0) {
            // Retrieve options from GUI:

            try {
                UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
            } catch (ClassNotFoundException | InstantiationException | UnsupportedLookAndFeelException | IllegalAccessException e) {
                Log.warning.println("Error setting cross-platform look and feel.");
            }

            try {
                SwingUtilities.invokeAndWait(() -> {
                    if (!getOptionsGUI(options))
                        System.exit(0);

                    setupGUIOutput();
                });
            } catch (InterruptedException | InvocationTargetException e) {
                e.printStackTrace();
            }


        } else {
            getCLIOptions(args, options);
        }

        // Run ACGAnnotator
        try {
            new ReMapTool(options);

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
