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
import org.apache.commons.math3.analysis.solvers.BrentSolver;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NoBracketingException;

import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import pitchfork.Pitchforks;

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
	public Input<RealParameter> popSizeAboveOriginInput = new Input<RealParameter>("popSizeAboveOrigin",
			"Constant above origin effective population sizes.", (RealParameter)null);
	public Input<String> fileNameInput = new Input<>("fileName",
			"Name of file to which gene trees will be written. Only if not using logger.");

	public Input<RealParameter> bottleneckStrengthInput = new Input<RealParameter>("bottleneckStrength",
			"Strength of the bottleneck in scaled time", Validate.REQUIRED);

	public Input<RealParameter> maxBoundInput = new Input<RealParameter>("maxBound",
			"Maximum value of time interval when computing next unobserved transmission time. Default 1000",
			Input.Validate.OPTIONAL);

	public Input<RealParameter> originInput = new Input<RealParameter>("origin", "origin",
			Input.Validate.REQUIRED);
	public Input<RealParameter> birthRateInput = new Input<RealParameter>("birthRate", "Birth rate",
			Input.Validate.REQUIRED);
	public Input<Function> deathRateInput = new Input<Function>("deathRate", "Death rate", Input.Validate.REQUIRED);
	public Input<RealParameter> samplingRateInput = new Input<RealParameter>("samplingRate",
			"Sampling rate per individual", Input.Validate.REQUIRED);

	public Input<RealParameter> samplingExtantRateInput = new Input<RealParameter>("samplingExtantRate",
			"Sampling rate of extant tips");

	public Input<Boolean> hiddenKnownInput = new Input<>("complete",
			"Simulate hidden bottlenecks", false);


	public Tree transmissionTree;
    public TraitSet sampleCounts;
	private int transmissionNodeCount;

	private double bottleneckStrength;
	private double maxBound;

	public int hiddenNodes;
	public List<Double> hiddenNodeTimes;

	private boolean hiddenKnown;
	public HashMap<String, Node> geneTreeEventAssignment = new HashMap<>();
	public HashMap<String, Node> geneTreeSampleAssignment;

	private double samplingExtantRate;

//	@Override
//	public void init(PrintStream out) {
//	}

	@Override
	public void close(PrintStream out) {
	}


    @Override
    public void initAndValidate() {

		transmissionTree = transmissionTreeInput.get();
		transmissionNodeCount = transmissionTree.getNodeCount();
		popSizesInput.get().setDimension(transmissionNodeCount);
		bottleneckStrength = bottleneckStrengthInput.get().getArrayValue(0);
		maxBound = maxBoundInput.get() != null ? maxBoundInput.get().getArrayValue(0) : 10000;

		samplingExtantRate = samplingExtantRateInput.get() == null ? 1 : samplingExtantRateInput.get().getArrayValue(0);

        sampleCounts = sampleCountsInput.get();

		hiddenNodes = 0;
		hiddenNodeTimes = new ArrayList();
		hiddenKnown = false;

		for (Node n : transmissionTree.getNodesAsArray()) {
			if (n.isLeaf() && sampleCounts.getValue(n.getID()) == 0) {
				hiddenKnown = true;
				if (n.getParent().getChild(0) != n)
					System.exit(0);
//				break;
			}
			if (hiddenKnownInput.get()) {
				hiddenKnown = true;
			}
		}

		Tree tree;
		while (true) {
			tree = getSimulatedGeneTree();
			List<Node> trueNodes = Pitchforks.getTrueInternalNodes(tree);
//			if (trueNodes.size() < 3 || (trueNodes.get(0).getHeight() != trueNodes.get(1).getHeight()
//					&& trueNodes.get(0).getHeight() != trueNodes.get(2).getHeight()
//					&& trueNodes.get(1).getHeight() != trueNodes.get(2).getHeight())) {
//				trueNodes = trueNodes.stream().sorted(Comparator.comparing(Node::getHeight))
//						.collect(Collectors.toList());
//				if (trueNodes.get(0).getHeight() == trueNodes.get(1).getHeight()
//						&& Pitchforks.getLogicalChildren(trueNodes.get(0)).size() == 2
//						&& Pitchforks.getLogicalChildren(trueNodes.get(1)).size() == 3) {
					assignFromWithoutID(tree);
					break;
//				}
//			}

//			}

//			if (trueNodes.size() == 3) {
//				trueNodes = trueNodes.stream().sorted(Comparator.comparing(Node::getHeight))
//						.collect(Collectors.toList());
//				List<Node> children = Pitchforks.getLogicalChildren(trueNodes.get(2));
//				if (children.size() == 2 && (children.get(0).getNr() == 4 || children.get(1).getNr() == 4)
//						&& Pitchforks.getLogicalChildren(trueNodes.get(0)).size() == 3) {
//					assignFromWithoutID(tree);
//					break;
//				}
//			}


//			if (trueNodes.size() == 2 && (Pitchforks.getLogicalChildren(trueNodes.get(0)).size() == 4
//					|| Pitchforks.getLogicalChildren(trueNodes.get(0)).size() == 5)) {
//				assignFromWithoutID(tree);
//				break;
//			}
//			else
//				tree = new Tree();
			
//			if (trueNodes.size() == 4) {
//				trueNodes = trueNodes.stream().sorted(Comparator.comparing(Node::getHeight))
//						.collect(Collectors.toList());
////			if (trueNodes.get(0).getHeight() == trueNodes.get(1).getHeight()
////						&& trueNodes.get(0).getHeight() != trueNodes.get(2).getHeight()
////						&& trueNodes.get(2).getHeight() != trueNodes.get(3).getHeight()
////						&& trueNodes.get(3).isRoot()) {
////				assignFromWithoutID(tree);
////				break;
////				}
//				if (trueNodes.get(0).getHeight() == trueNodes.get(1).getHeight()
//						&& trueNodes.get(0).getHeight() != trueNodes.get(2).getHeight()
//						&& trueNodes.get(2).getHeight() != trueNodes.get(3).getHeight()
//						&& trueNodes.get(3).isRoot()
//						&& (trueNodes.get(3).getChild(0).getNr() == 5
//								| trueNodes.get(3).getChild(1).getNr() == 5)
//						&& (Pitchforks.getGroup(trueNodes.get(0)).size() == 1
//								| Pitchforks.getGroup(trueNodes.get(1)).size() == 1)) {
//				assignFromWithoutID(tree);
//				break;
//				}
//			}

//			tree = new Tree();

		}

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

		geneTreeSampleAssignment = new HashMap<>();

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
		Node aboveOrigin = new Node();
		if (popSizeAboveOriginInput.get()==null)
			popSize.put(aboveOrigin, 1.0);
		else
			popSize.put(aboveOrigin, popSizeAboveOriginInput.get().getValue(0));

		while (getTotalLineageCount(activeLineages) > 1 || !sortedTransmissionTreeNodes.isEmpty()) {

			Node minCoal = new Node();
			Node minNonObsTr = new Node();
			double minCoalTime = Double.POSITIVE_INFINITY;
			double minNonObsTrTime = Double.POSITIVE_INFINITY;

			for (Node transmissionNode : activeLineages.keySet()) {
				if ((!transmissionNode.isRoot() && transmissionNode.getParent().getHeight() > t)
						|| transmissionNode.isRoot()) {
					int k = activeLineages.get(transmissionNode).size();
					double Ne = popSize.get(transmissionNode);
					if (transmissionNode.isRoot()) {
						if (k > 1) {
							double NeAbove = popSize.get(aboveOrigin);
							double timeToNextCoalOrig = Randomizer
									.nextExponential(0.5 * k * (k - 1) * 1.0 / NeAbove);

							if (t < originInput.get().getArrayValue(0))
								timeToNextCoalOrig += originInput.get().getArrayValue(0) - t;
							double timeToNextCoalBefore = Randomizer
									.nextExponential(0.5 * k * (k - 1) * 1.0 / Ne);
							if (t + timeToNextCoalBefore > originInput.get().getArrayValue(0))
								timeToNextCoalBefore = Double.POSITIVE_INFINITY;

							double timeNextCoalMin = Math.min(timeToNextCoalOrig, timeToNextCoalBefore);

							if (timeNextCoalMin < minCoalTime) {
								minCoalTime = timeNextCoalMin;
								minCoal = transmissionNode;
							}

						}
					} else {
						if (t > originInput.get().getArrayValue(0)) {
							Ne = popSize.get(aboveOrigin);
						}
						if (k > 1) {
							double timeToNextCoal = Randomizer
									.nextExponential(0.5 * k * (k - 1) * 1.0 / Ne);
							if (timeToNextCoal < minCoalTime) {
								minCoalTime = timeToNextCoal;
								minCoal = transmissionNode;
							}

						}
					}



					if (!hiddenKnown) {
						double endTime = transmissionNode.isRoot() ? maxBound
								: transmissionNode.getParent().getHeight();
//						for (int i = 0; i < 20; i++) {
//							System.out.println(soveForTime(t,
//									birthRateInput.get().getArrayValue(0), deathRateInput.get().getArrayValue(0),
//									samplingRateInput.get().getArrayValue(0), endTime) - t);
//						}
//						System.out.println("#####");
//						for (int i = 0; i < 20; i++) {
//							System.out.println(soveForTime(0,
//									birthRateInput.get().getArrayValue(0), deathRateInput.get().getArrayValue(0),
//									samplingRateInput.get().getArrayValue(0), endTime - t));
//						}
//						double timeToNextNonObsTr = soveForTime(0,
//								birthRateInput.get().getArrayValue(0), deathRateInput.get().getArrayValue(0),
//								samplingRateInput.get().getArrayValue(0), endTime - t);
						double timeToNextNonObsTr = soveForTime(t,
								birthRateInput.get().getArrayValue(0), deathRateInput.get().getArrayValue(0),
								samplingRateInput.get().getArrayValue(0), endTime) - t;
//						double timeToNextNonObsTr = Double.POSITIVE_INFINITY;
//						if (transmissionNode.isRoot())
//							timeToNextNonObsTr = Randomizer
//									.nextExponential(birthRateInput.get().getArrayValue(0) * integralP_0(
//											transmissionNode.getHeight(), originInput.get().getArrayValue(0),
//											birthRateInput.get().getArrayValue(0),
//											deathRateInput.get().getArrayValue(0),
//											samplingRateInput.get().getArrayValue(0)));
//						else
//							timeToNextNonObsTr = Randomizer
//									.nextExponential(birthRateInput.get().getArrayValue(0) * integralP_0(
//											transmissionNode.getHeight(), transmissionNode.getParent().getHeight(),
//											birthRateInput.get().getArrayValue(0),
//											deathRateInput.get().getArrayValue(0),
//											samplingRateInput.get().getArrayValue(0)));

						if (timeToNextNonObsTr + t < originInput.get().getArrayValue(0)
								&& timeToNextNonObsTr < minNonObsTrTime) {
							minNonObsTrTime = timeToNextNonObsTr;
							minNonObsTr = transmissionNode;
						}
					}
				}
			}

			
			double dt = Math.min(minCoalTime, minNonObsTrTime);

			if (!sortedTransmissionTreeNodes.isEmpty() && t + dt >= sortedTransmissionTreeNodes.get(0).getHeight()) {
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
						geneTreeEventAssignment.put(geneTreeSampleNode.getID(), transmissionNode);
						geneTreeSampleAssignment.put(geneTreeSampleNode.getID(), transmissionNode);
                    }

                } else {
					// Observed transmission

					Node recipient = transmissionNode.getChildren().get(1);
					Node donor = transmissionNode.getChild(0);
					if (donor.isLeaf() && sampleCounts.getValue(donor.getID()) == 0) {
						hiddenNodes += 1;
						hiddenNodeTimes.add(donor.getHeight());
					}
					List<Node> recipientLineages = activeLineages.get(recipient);
					double recipientNe = popSize.get(recipient);

					// From Tim: Can this be extracted to a method?  More-or-less duplicated below.
					if  (bottleneckStrength > 0 && recipientLineages.size() > 1 &&
								 !recipient.isDirectAncestor()) {
						double duplicateTime = t;
						double stopTime = t + bottleneckStrength;
						while (duplicateTime < stopTime) {
							int nrLineages = recipientLineages.size();
							double deltaT = Randomizer
									.nextExponential(nrLineages * (nrLineages - 1) * 0.5 * 1.0 / recipientNe);
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
								
								geneTreeEventAssignment.put(parent.getID(), recipient);
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
					
					geneTreeEventAssignment.put(parent.getID(), minCoal);

				} else if (minNonObsTrTime == dt) {
					// Unobserved transmission
					hiddenNodes += 1;
					hiddenNodeTimes.add(minNonObsTrTime);

					List<Node> lineageList = activeLineages.get(minNonObsTr);
					double Ne = popSize.get(minNonObsTr);

					// Duplicate of above code
					if (bottleneckStrength > 0 && lineageList.size() > 1) {
						double duplicateTime = t;
						double stopTime = t + bottleneckStrength;
						while (duplicateTime < stopTime) {
							int nrLineages = lineageList.size();
							double deltaT = Randomizer
									.nextExponential(nrLineages * (nrLineages - 1) * 0.5 * 1.0 / Ne);
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

								geneTreeEventAssignment.put(parent.getID(), minNonObsTr);
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
		Tree tree = new Tree(activeLineages.get(transmissionTree.getRoot()).get(0));

		return tree;
//		String newick = "(((((((t7_19:0.07166212308424236,t7_14:0.07166212308424236)170:0.14445426969305297,t7_114:0.21611639277729533)189:0.2956333653688481,((t7_17:0.027740724119904574,t7_15:0.027740724119904574)144:0.4593348963341476,((t7_113:0.07264985101614894,(t7_116:0.049326480975383495,(t7_118:0.00466680201469649,t7_111:0.00466680201469649)126:0.044659678960687)152:0.02332337004076545)171:0.2330343700854443,(((t7_110:0.04734176224658845,(t7_119:0.008674367820234686,t7_13:0.008674367820234686)130:0.03866739442635376)150:0.008880767187851739,(t7_112:0.023358392466842376,t7_11:0.023358392466842376)140:0.03286413696759781)167:0.07179055243214084,t7_117:0.12801308186658103)181:0.17767113923501224)194:0.1813913993524589)207:0.024674137692091258)209:0.09857747089453173,((t7_120:0.07520889922830291,t7_18:0.07520889922830291)172:0.10965799054716378,((t7_115:0.01601751226081743,t7_12:0.01601751226081743)135:0.022579766639740314,t7_16:0.038597278900557745)148:0.14626961087490895)185:0.42546033926520843)211:0.4535720669754433,(((((t2_13:0.02918484707520955,t2_116:0.02918484707520955)145:0.020141633900173945,t2_18:0.049326480975383495)153:0.0,t2_111:0.049326480975383495)154:0.0,(t2_110:3.419806080766728E-5,t2_15:3.419806080766728E-5)120:0.049292282914575825)157:0.0,((t2_118:0.023510000768121597,t2_113:0.023510000768121597)141:0.025816480207261898,t2_16:0.049326480975383495)158:0.0)160:1.014572815040735)215:0.4670127663389929,(((((((t2_115:0.007176066687872195,t2_14:0.007176066687872195)127:0.0084520935416071,t2_117:0.015628160229479295)134:0.033577076922281596,t2_11:0.049205237151760894)151:1.2124382362260072E-4,(t2_112:0.007808194975578015,t2_120:0.007808194975578015)129:0.041518285999805477)155:0.0,t2_12:0.049326480975383495)156:0.0,(t2_114:6.523965399714786E-4,t2_19:6.523965399714786E-4)122:0.04867408443541202)159:0.0,(t2_17:0.049326480975383495,t2_119:0.049326480975383495)161:0.0)162:1.4815855813797278)233:0.24456990037562476,((((t22_115:0.07692848329536672,(t22_12:0.04826563114346327,t22_11:0.04826563114346327)220:0.02866285215190345)223:0.08774574008069558,(t22_15:0.03279740493763805,(t22_119:0.001604184487313809,t22_19:0.001604184487313809)216:0.031193220450324244)217:0.13187681843842425)228:0.34003632968403963,(t22_113:0.5047105530601019,t22_116:0.5047105530601019)235:0.0)236:0.0,((t22_120:0.15946155320275102,(((t6_116:0.12351858292000446,(t6_118:0.06682494716315725,t6_15:0.06682494716315725)169:0.056693635756847216)180:1.1375525238978321,(((((((t11_16:0.23685842861326617,t11_11:0.23685842861326617)190:0.17585008388570267,((t11_14:0.36281483805464576,t11_120:0.36281483805464576)202:0.04622730203211073,((t11_112:0.025058232964184402,t11_19:0.025058232964184402)143:0.08993147663358961,t11_119:0.114989709597774)179:0.2940524304889825)203:0.0036663724122123487)204:0.029190267692703886,((((t11_114:0.020617892442161778,t11_110:0.020617892442161778)137:0.002418086887048008,t11_13:0.023035979329209786)139:0.1717929017585616,(t11_116:0.02156298961542966,t11_17:0.02156298961542966)138:0.17326589147234173)186:0.05686643395816046,(t11_15:0.009179012235472336,t11_18:0.009179012235472336)131:0.2425163028104595)192:0.19020346514574088)205:0.059500190183397916,(t6_110:0.13562229534779563,((t6_113:0.004393866360303908,t6_14:0.004393866360303908)125:0.12747836932175716,t6_12:0.13187223568206105)182:0.0037500596657345786)183:0.36577667502727634)208:0.09072397989883974,(((t6_13:0.10138171547623469,t6_18:0.10138171547623469)177:0.15479254686282856,t6_19:0.25617426233906326)193:0.07937240023936531,(t6_111:0.10649369622985964,t6_120:0.10649369622985964)178:0.22905296634856892)201:0.25657628769548313)210:0.018748872135805006,(t6_17:0.24253649440064937,(t6_114:0.089624442597338,(t6_119:0.017940731505720117,(t6_112:5.595281270931729E-4,t6_115:5.595281270931729E-4)121:0.017381203378626944)136:0.07168371109161788)174:0.15291205180331136)191:0.3683353280090673)212:0.2855188111645959,(((((t11_118:0.007714609585906115,t11_111:0.007714609585906115)128:0.04759744378162702,t11_117:0.055312053367533136)165:0.1556004877564606,((t11_12:0.0012074927478597603,t11_115:0.0012074927478597603)123:0.0227065922190829,t11_113:0.02391408496694266)142:0.18699845615705107)188:0.23098623906767898,((t10_13:0.32141527771213685,(((t10_11:0.03786853405480256,t10_120:0.03786853405480256)147:0.27394918796457335,(t10_14:0.17820702916117653,(((t10_18:0.014372408453372577,t10_118:0.014372408453372577)133:0.041834573677790164,t10_110:0.05620698213116274)166:0.04386314755917801,((t10_19:0.0018458490174113984,t10_16:0.0018458490174113984)124:0.09507565783212625,t10_116:0.09692150684953765)175:0.0031486228408030975)176:0.07813689947083578)184:0.13361069285819938)196:0.00959755569276094,(t10_115:0.057901488462270605,(t10_17:0.05520365717497741,t10_113:0.05520365717497741)164:0.002697831287293194)168:0.26351378924986624)198:0.0)199:0.0,(((t10_119:0.029608432204647848,t10_117:0.029608432204647848)146:0.28181890999952774,(t10_111:0.19809026948150044,t10_15:0.19809026948150044)187:0.11333707272267518)195:0.009987935507961232,((t10_112:0.009622801285743597,t10_12:0.009622801285743597)132:0.040021840342847274,t10_114:0.04964464162859087)163:0.27177063608354596)197:0.0)200:0.12048350247953588)206:0.21556477592823176,((t6_16:0.04051307477190871,t6_117:0.04051307477190871)149:0.0418951516526287,t6_11:0.0824082264245374)173:0.5750553296953684)213:0.23892707745440678)214:0.36468047324352404)222:0.03570772139389078,t22_118:0.10719471506912548)224:0.05226683813362554)227:0.2585281281995435,(((((t22_111:0.1137384196592508,(((t22_117:0.0416246318055451,t22_112:0.0416246318055451)218:0.0036839569488305024,t22_14:0.0453085887543756)219:0.014153584091269433,t22_17:0.059462172845645034)221:0.054276246813605766)225:0.023648796066987998,t22_16:0.1373872157262388)226:0.029951586546849596,t22_18:0.1673388022730884)229:0.0047256445098207145,t22_114:0.1720644467829091)230:0.09117028498199975,(t22_110:0.23761604307637008,t22_13:0.23761604307637008)231:0.025618688688538782)232:0.15475494963738567)234:0.0867208716578074)237:0.08118729652802825)238:0.0";
//		TreeParser tr = new TreeParser();
//		tr.initByName("newick", newick, "adjustTipHeights", false, "IsLabelledNewick", true);
//		return tr;
    }

	// Solve an inverse transform problem for time t of next occurence
	private double soveForTime(double currentTime, double lambda, double mu, double psi, double endTime) {

		double u = Randomizer.nextDouble();
		double c_1 = Math.abs(Math.sqrt( Math.pow((lambda - mu - psi), 2) + 4*lambda*psi ));
		double c_2 = -(lambda - mu - 2 * lambda * samplingExtantRate - psi) / c_1;
		
//		NewtonRaphsonSolver solver = new NewtonRaphsonSolver(1E-10);

		UnivariateDifferentiableFunction f = new UnivariateDifferentiableFunction() {

			private double function(double t) {
				return ((-2 * Math.log(Math.abs((c_2 - 1) * Math.exp(-c_1 * t) - c_2 - 1))
						+ (-c_1 * t) + (mu + psi + lambda) * t) / 2);

			}

			private DerivativeStructure ds(DerivativeStructure t) {
				return ((t.multiply(-c_1)).exp().multiply(c_2 - 1).subtract(c_2 + 1).abs().log().multiply(-2)
						.add(t.multiply(-c_1).exp().log()).add(t.multiply(lambda + mu + psi))).divide(2);
			}

			// solve fot t: int_{currentTime}^{t} f(x) = -ln(1-u)
			@Override
			public double value(double t) {
				return this.function(t) + Math.log(1-u) - this.function(currentTime);
			}

			@Override
			public DerivativeStructure value(DerivativeStructure t) throws DimensionMismatchException {

				return this.ds(t).add(Math.log(1 - u)).subtract(this.ds(new DerivativeStructure(1, 1, currentTime)));

			}
		};

		BrentSolver solver = new BrentSolver(1e-15);
		double time = Double.POSITIVE_INFINITY; // return infinity if there are no roots within time interval of the
												// transmission tree
												// branch

		double kls = 0.0;
		double klss = Math.log(1 - u);
		try {
			time = solver.solve(100000, f, currentTime, endTime);
		} catch (NoBracketingException e) {
		}
		if (time != Double.POSITIVE_INFINITY && (time >= endTime || time <= currentTime))
			time = Double.POSITIVE_INFINITY;
		if (time != Double.POSITIVE_INFINITY)
			kls = integralP_0(currentTime, time, lambda, mu, psi);

		return time;
		
//		return solver.solve(100000, f, currentTime, maxBound);

	}

	private double integralP_0(double t_0, double t_1, double lambda, double mu, double psi) {
		double c1 = Math.abs(Math.sqrt(Math.pow((lambda - mu - psi), 2) + 4 * lambda * psi));
		double c2 = -(lambda - mu - 2 * lambda * samplingExtantRate - psi) / c1;

		double ans = (1.0 / (2.0 * lambda)) * ((t_1 - t_0) * (mu + psi + lambda - c1) + 2.0 * Math
				.log(((c2 - 1) * Math.exp(-c1 * t_0) - c2 - 1) / ((c2 - 1) * Math.exp(-c1 * t_1) - c2 - 1)));

		return ans;
	}

}
