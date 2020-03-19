package tnt.likelihood;



import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import tnt.simulator.SimulatedGeneTree;

/**
 *
 * @author Ugne Stolz
 * 
 */
@Description("Extracts the intervals from a gene tree. Points in the intervals " +
        "are defined by the heights of nodes in the tree.")
public class GeneTreeIntervals extends CalculationNode {
	public Input<Tree> geneTreeInput = new Input<Tree>("geneTree", "Gene tree for which to calculate the intervals");

	public Input<SimulatedGeneTree> simulatedGeneTreeInput=new Input<SimulatedGeneTree>("simulatedGeneTree","Gene tree for which to calculate the intervals",
			Validate.XOR, geneTreeInput);

	public Input<Tree> transmissionTreeInput = new Input<>("transmissionTreeInput",
			"Fully labeled transmission tree on which to simulate gene trees",
			Input.Validate.REQUIRED);

	private SimulatedGeneTree geneTree;
    
	private List<GeneTreeEvent> geneTreeEventList, storedGeneTreeEventList;
	HashMap<Node, Integer> activeLineagesPerTransmissionTreeNode;
	HashMap<Node, List<GeneTreeEvent>> eventsPerTransmissionTreeNode;

	HashMap<Node, List<String>> nodesPerTransmissionTreeNode;

	HashMap<String, Node> geneTreeNodeAssignment;
	public boolean eventListDirty = true;

	@Override
	public void initAndValidate() {
		geneTree = simulatedGeneTreeInput.get();

//		String newick = "(((((((t7_19:0.07166212308424236,t7_14:0.07166212308424236)170:0.14445426969305297,t7_114:0.21611639277729533)189:0.2956333653688481,((t7_17:0.027740724119904574,t7_15:0.027740724119904574)144:0.4593348963341476,((t7_113:0.07264985101614894,(t7_116:0.049326480975383495,(t7_118:0.00466680201469649,t7_111:0.00466680201469649)126:0.044659678960687)152:0.02332337004076545)171:0.2330343700854443,(((t7_110:0.04734176224658845,(t7_119:0.008674367820234686,t7_13:0.008674367820234686)130:0.03866739442635376)150:0.008880767187851739,(t7_112:0.023358392466842376,t7_11:0.023358392466842376)140:0.03286413696759781)167:0.07179055243214084,t7_117:0.12801308186658103)181:0.17767113923501224)194:0.1813913993524589)207:0.024674137692091258)209:0.09857747089453173,((t7_120:0.07520889922830291,t7_18:0.07520889922830291)172:0.10965799054716378,((t7_115:0.01601751226081743,t7_12:0.01601751226081743)135:0.022579766639740314,t7_16:0.038597278900557745)148:0.14626961087490895)185:0.42546033926520843)211:0.4535720669754433,(((((t2_13:0.02918484707520955,t2_116:0.02918484707520955)145:0.020141633900173945,t2_18:0.049326480975383495)153:0.0,t2_111:0.049326480975383495)154:0.0,(t2_110:3.419806080766728E-5,t2_15:3.419806080766728E-5)120:0.049292282914575825)157:0.0,((t2_118:0.023510000768121597,t2_113:0.023510000768121597)141:0.025816480207261898,t2_16:0.049326480975383495)158:0.0)160:1.014572815040735)215:0.4670127663389929,(((((((t2_115:0.007176066687872195,t2_14:0.007176066687872195)127:0.0084520935416071,t2_117:0.015628160229479295)134:0.033577076922281596,t2_11:0.049205237151760894)151:1.2124382362260072E-4,(t2_112:0.007808194975578015,t2_120:0.007808194975578015)129:0.041518285999805477)155:0.0,t2_12:0.049326480975383495)156:0.0,(t2_114:6.523965399714786E-4,t2_19:6.523965399714786E-4)122:0.04867408443541202)159:0.0,(t2_17:0.049326480975383495,t2_119:0.049326480975383495)161:0.0)162:1.4815855813797278)233:0.24456990037562476,((((t22_115:0.07692848329536672,(t22_12:0.04826563114346327,t22_11:0.04826563114346327)220:0.02866285215190345)223:0.08774574008069558,(t22_15:0.03279740493763805,(t22_119:0.001604184487313809,t22_19:0.001604184487313809)216:0.031193220450324244)217:0.13187681843842425)228:0.34003632968403963,(t22_113:0.5047105530601019,t22_116:0.5047105530601019)235:0.0)236:0.0,((t22_120:0.15946155320275102,(((t6_116:0.12351858292000446,(t6_118:0.06682494716315725,t6_15:0.06682494716315725)169:0.056693635756847216)180:1.1375525238978321,(((((((t11_16:0.23685842861326617,t11_11:0.23685842861326617)190:0.17585008388570267,((t11_14:0.36281483805464576,t11_120:0.36281483805464576)202:0.04622730203211073,((t11_112:0.025058232964184402,t11_19:0.025058232964184402)143:0.08993147663358961,t11_119:0.114989709597774)179:0.2940524304889825)203:0.0036663724122123487)204:0.029190267692703886,((((t11_114:0.020617892442161778,t11_110:0.020617892442161778)137:0.002418086887048008,t11_13:0.023035979329209786)139:0.1717929017585616,(t11_116:0.02156298961542966,t11_17:0.02156298961542966)138:0.17326589147234173)186:0.05686643395816046,(t11_15:0.009179012235472336,t11_18:0.009179012235472336)131:0.2425163028104595)192:0.19020346514574088)205:0.059500190183397916,(t6_110:0.13562229534779563,((t6_113:0.004393866360303908,t6_14:0.004393866360303908)125:0.12747836932175716,t6_12:0.13187223568206105)182:0.0037500596657345786)183:0.36577667502727634)208:0.09072397989883974,(((t6_13:0.10138171547623469,t6_18:0.10138171547623469)177:0.15479254686282856,t6_19:0.25617426233906326)193:0.07937240023936531,(t6_111:0.10649369622985964,t6_120:0.10649369622985964)178:0.22905296634856892)201:0.25657628769548313)210:0.018748872135805006,(t6_17:0.24253649440064937,(t6_114:0.089624442597338,(t6_119:0.017940731505720117,(t6_112:5.595281270931729E-4,t6_115:5.595281270931729E-4)121:0.017381203378626944)136:0.07168371109161788)174:0.15291205180331136)191:0.3683353280090673)212:0.2855188111645959,(((((t11_118:0.007714609585906115,t11_111:0.007714609585906115)128:0.04759744378162702,t11_117:0.055312053367533136)165:0.1556004877564606,((t11_12:0.0012074927478597603,t11_115:0.0012074927478597603)123:0.0227065922190829,t11_113:0.02391408496694266)142:0.18699845615705107)188:0.23098623906767898,((t10_13:0.32141527771213685,(((t10_11:0.03786853405480256,t10_120:0.03786853405480256)147:0.27394918796457335,(t10_14:0.17820702916117653,(((t10_18:0.014372408453372577,t10_118:0.014372408453372577)133:0.041834573677790164,t10_110:0.05620698213116274)166:0.04386314755917801,((t10_19:0.0018458490174113984,t10_16:0.0018458490174113984)124:0.09507565783212625,t10_116:0.09692150684953765)175:0.0031486228408030975)176:0.07813689947083578)184:0.13361069285819938)196:0.00959755569276094,(t10_115:0.057901488462270605,(t10_17:0.05520365717497741,t10_113:0.05520365717497741)164:0.002697831287293194)168:0.26351378924986624)198:0.0)199:0.0,(((t10_119:0.029608432204647848,t10_117:0.029608432204647848)146:0.28181890999952774,(t10_111:0.19809026948150044,t10_15:0.19809026948150044)187:0.11333707272267518)195:0.009987935507961232,((t10_112:0.009622801285743597,t10_12:0.009622801285743597)132:0.040021840342847274,t10_114:0.04964464162859087)163:0.27177063608354596)197:0.0)200:0.12048350247953588)206:0.21556477592823176,((t6_16:0.04051307477190871,t6_117:0.04051307477190871)149:0.0418951516526287,t6_11:0.0824082264245374)173:0.5750553296953684)213:0.23892707745440678)214:0.36468047324352404)222:0.03570772139389078,t22_118:0.10719471506912548)224:0.05226683813362554)227:0.2585281281995435,(((((t22_111:0.1137384196592508,(((t22_117:0.0416246318055451,t22_112:0.0416246318055451)218:0.0036839569488305024,t22_14:0.0453085887543756)219:0.014153584091269433,t22_17:0.059462172845645034)221:0.054276246813605766)225:0.023648796066987998,t22_16:0.1373872157262388)226:0.029951586546849596,t22_18:0.1673388022730884)229:0.0047256445098207145,t22_114:0.1720644467829091)230:0.09117028498199975,(t22_110:0.23761604307637008,t22_13:0.23761604307637008)231:0.025618688688538782)232:0.15475494963738567)234:0.0867208716578074)237:0.08118729652802825)238:0.0";
//		TreeParser tr = new TreeParser();
//		tr.initByName("newick", newick, "adjustTipHeights", false, "IsLabelledNewick", true);
//		SimulatedGeneTree trr = new SimulatedGeneTree();

//		geneTree = tr.copy();

		storedGeneTreeEventList = new ArrayList<>();
		activeLineagesPerTransmissionTreeNode = new HashMap<>();
		eventsPerTransmissionTreeNode = new HashMap<>();

	}

	HashMap<Node, List<GeneTreeEvent>> getGeneTreeEventList() {
		update();
		return eventsPerTransmissionTreeNode;
	}

	private void update() {
        if (!eventListDirty)
            return;

		geneTreeNodeAssignment = geneTree.geneTreeSampleAssignment;
		HashMap<Double, List<Node>> fakeBifurcations = new HashMap<>();
		HashMap<Node, List<Node>> fakeBifurcations2 = new HashMap<>();
		HashMap<Double, List<Node>> nodeTime = new HashMap<>();
		geneTreeEventList = new ArrayList<GeneTreeEvent>();
		activeLineagesPerTransmissionTreeNode = new HashMap<>();
		eventsPerTransmissionTreeNode = new HashMap<>();

//		int id = 0;
//		for (Node trNode : transmissionTreeInput.get().getNodesAsArray()) {
//			if (!trNode.isLeaf()) {
//				Node mockNode = new Node();
//				mockNode.setID("mock_" + id);
//				id++;
//				mockNode.setHeight(trNode.getHeight());
//				GeneTreeEvent startEvent = new GeneTreeEvent();
//				startEvent.node = mockNode;
//				startEvent.type = GeneTreeEvent.GeneTreeEventType.TANSMISSION;
//				startEvent.time = trNode.getHeight();
//				geneTreeEventList.add(startEvent);
//				geneTreeNodeAssignment.put(startEvent.node.getID(), trNode);
////				geneTree.geneTreeEventAssignment.put(startEvent.node.getID(), trNode);
//			}
//		}

		for (Node n : geneTree.getNodesAsArray()) {

			if (!n.isRoot() && n.getParent().getHeight() - n.getHeight() == 0) {
				Node multifurcationParent = getMultifurcationParent(n, n.getParent());
				if (fakeBifurcations2.containsKey(multifurcationParent))
					fakeBifurcations2.get(multifurcationParent).add(n);
				else {
					List<Node> list = new ArrayList<Node>();
					list.add(n);
					fakeBifurcations2.put(multifurcationParent, list);
				}
			}

			Node trNode = null;
			if (!n.isLeaf()) {
				Node child = n.getChild(0);
				trNode = geneTreeNodeAssignment.get(child.getID());
				while (trNode == null) {
					child = child.getChild(0);
					trNode = geneTreeNodeAssignment.get(child.getID());
				}
				while (!trNode.isRoot() && n.getHeight() > trNode.getParent().getHeight()) {
					trNode = trNode.getParent();
				}
				geneTreeNodeAssignment.put(n.getID(), trNode);
			}
		}

		nodeTime = new HashMap<Double, List<Node>>();
		for (Node node : geneTree.getNodesAsArray()) {
			if (node.isRoot() || (!node.isLeaf() && node.getParent().getHeight() != node.getHeight())) {
				if (nodeTime.containsKey(node.getHeight())) {
					nodeTime.get(node.getHeight()).add(node);
				}
				else {
					List<Node> list = new ArrayList<Node>();
					list.add(node);
					nodeTime.put(node.getHeight(), list);
				}
			}
		}
		
		for (Double time : nodeTime.keySet()) {
			Node first = nodeTime.get(time).get(0);

			GeneTreeEvent event = new GeneTreeEvent();
			event.time = first.getHeight();
			event.node = first;
			event.fakeBifCount = 0;
			if (fakeBifurcations2.containsKey(first)) {
				event.fakeBifCount += fakeBifurcations2.get(first).size();
				event.multiCoalCount += 1;
				event.multiCoalSize.add(fakeBifurcations2.get(first).size() + 2);
				event.type = GeneTreeEvent.GeneTreeEventType.MULTIFURCATION;
//			}

			if (nodeTime.get(time).size() > 1) {

				for (Node rest : nodeTime.get(time)) {
					if (first != rest &&
							geneTreeNodeAssignment.get(first.getID()) == geneTreeNodeAssignment.get(rest.getID())) {
							if (fakeBifurcations2.containsKey(rest)) {
								event.fakeBifCount += fakeBifurcations2.get(rest).size();
								event.multiCoalSize.add(fakeBifurcations2.get(rest).size() + 2);
							} else {
								event.multiCoalSize.add(2);
							}
							event.multiCoalCount += 1;
						event.type = GeneTreeEvent.GeneTreeEventType.MULTIFURCATION;
					}
				}
				}
			} else if (first.getChildCount() == 2)
				event.type = GeneTreeEvent.GeneTreeEventType.BIFURCATION;

			geneTreeEventList.add(event);
		}

		for (Node n : geneTree.getNodesAsArray()) {
			if (n.getChildCount() == 0) {
				GeneTreeEvent event = new GeneTreeEvent();
				event.time = n.getHeight();
				event.node = n;
				event.type = GeneTreeEvent.GeneTreeEventType.SAMPLE;
				geneTreeEventList.add(event);
			}
		}




//		for (Node n : geneTree.getNodesAsArray()) {
//			if (!n.isRoot() && n.getParent().getHeight() - n.getHeight() == 0) {
//				if (fakeBifurcations.containsKey(n.getHeight()))
//					fakeBifurcations.get(n.getHeight()).add(n);
//				else {
//					List<Node> list = new ArrayList<Node>();
//					list.add(n);
//					fakeBifurcations.put(n.getHeight(), list);
//				}
//				continue;
//			}
//
//			GeneTreeEvent event = new GeneTreeEvent();
//			event.time = n.getHeight();
//			event.node = n;
//
//			if (n.getChildCount() == 0)
//				event.type = GeneTreeEvent.GeneTreeEventType.SAMPLE;
//			else if (n.getChildCount() == 2)
//				event.type = GeneTreeEvent.GeneTreeEventType.BIFURCATION;
//			else
//				throw new RuntimeException("Network node has illegal number of children.");
//
//			geneTreeEventList.add(event);
//		}
//		
//		HashMap<Double, List<GeneTreeEvent>> multiCoal = new HashMap<>();
//
//		for (GeneTreeEvent event : geneTreeEventList) {
//			if (fakeBifurcations.containsKey(event.time)) {
//				event.fakeBifCount = fakeBifurcations.get(event.time).size();
//				event.type = GeneTreeEvent.GeneTreeEventType.MULTIFURCATION;
//				event.multiCoalCount += 1;
//
//				if (multiCoal.containsKey(event.time))
//					multiCoal.get(event.time).add(event);
//				else {
//					List<GeneTreeEvent> coalList = new ArrayList<>();
//					coalList.add(event);
//					multiCoal.put(event.time, coalList);
//				}
//			}
//		}
//		
//		for (double time : multiCoal.keySet()) {
//			multiCoal.get(time).get(0).multiCoalSize.add(coalSize(multiCoal.get(time).get(0).node));
//			if (multiCoal.get(time).size() > 1) {
//				for (int i = 1; i < multiCoal.get(time).size(); i++) {
//					multiCoal.get(time).get(0).multiCoalCount += multiCoal.get(time).get(i).multiCoalCount;
//					multiCoal.get(time).get(0).multiCoalSize.add(coalSize(multiCoal.get(time).get(i).node));
//					geneTreeEventList.remove(multiCoal.get(time).get(i));
//				}
//			}
//		}



		geneTreeEventList = geneTreeEventList.stream()
				.sorted(Comparator.comparingDouble(e -> e.time))
				.collect(Collectors.toList());

//		System.out.println(this.geneTree.getRoot().toNewick());
        for (GeneTreeEvent event : geneTreeEventList) {

//			Node trNode = null;
//        	if (!event.node.isLeaf()) {
//				Node child = event.node.getChild(0);
//				trNode = geneTreeNodeAssignment.get(child.getID());
//				while (trNode == null) {
//					child = child.getChild(0);
//					trNode = geneTreeNodeAssignment.get(child.getID());
//				}
//				while (!trNode.isRoot() && event.node.getHeight() > trNode.getParent().getHeight()) {
//					trNode = trNode.getParent();
//				}
//
//				geneTreeNodeAssignment.put(event.node.getID(), trNode);
//			} else {
			Node trNode = geneTreeNodeAssignment.get(event.node.getID());
//			}

//			Node trNode = geneTree.geneTreeEventAssignment.get(event.node.getID());
			if (trNode.getID() == null || !trNode.getID().equals("t7_1"))
				continue;
			int nrLineage = activeLineagesPerTransmissionTreeNode.get(trNode) == null ? 0
					: activeLineagesPerTransmissionTreeNode.get(trNode);
			if (!trNode.isLeaf() && nrLineage == 0) {
				activeLineagesPerTransmissionTreeNode.put(trNode, getLineagesRecurse(trNode));
				nrLineage = activeLineagesPerTransmissionTreeNode.get(trNode);
			}

        	switch(event.type) {
			case SAMPLE:
				nrLineage += 1;
				activeLineagesPerTransmissionTreeNode.put(trNode, nrLineage);
				if (trNode.getParent().isFake())
					activeLineagesPerTransmissionTreeNode.put(trNode.getParent(), nrLineage);

				break;
			case BIFURCATION:
				nrLineage -= 1;
				activeLineagesPerTransmissionTreeNode.put(trNode, nrLineage);
//				if (!trNode.isRoot() && event.node.getHeight() == trNode.getParent().getHeight()) {
//					Node otherChild = trNode == trNode.getParent().getChild(0) ? trNode.getParent().getChild(1)
//							: trNode.getParent().getChild(0);
//					activeLineagesPerTransmissionTreeNode.put(trNode.getParent(),
//							nrLineage + activeLineagesPerTransmissionTreeNode.get(otherChild));
//				}
				break;
			case MULTIFURCATION:
				nrLineage -= (event.fakeBifCount + event.multiCoalCount);
				activeLineagesPerTransmissionTreeNode.put(trNode, nrLineage);
//				if (!trNode.isRoot() && event.node.getHeight() == trNode.getParent().getHeight()) {
//					Node otherChild = trNode == trNode.getParent().getChild(0) ? trNode.getParent().getChild(1)
//							: trNode.getParent().getChild(0);
//					activeLineagesPerTransmissionTreeNode.put(trNode.getParent(),
//							nrLineage + activeLineagesPerTransmissionTreeNode.get(otherChild));
//				}
				break;
        	}

			if (nrLineage <= 0)
				System.out.println();
			event.lineages = nrLineage;

			if (eventsPerTransmissionTreeNode.get(trNode) == null) {
				List<GeneTreeEvent> l = new ArrayList<>();
				l.add(event);
				eventsPerTransmissionTreeNode.put(trNode, l);
			} else {
				eventsPerTransmissionTreeNode.get(trNode).add(event);
			}
        }

		eventListDirty = false;
	}

	int getLineagesRecurse(Node trNode) {
		int nLineages = 0;
		if (activeLineagesPerTransmissionTreeNode.get(trNode.getChild(0)) == null)
			nLineages += getLineagesRecurse(trNode.getChild(0));
		else
			nLineages += activeLineagesPerTransmissionTreeNode.get(trNode.getChild(0));
		if (activeLineagesPerTransmissionTreeNode.get(trNode.getChild(1)) == null)
			nLineages += getLineagesRecurse(trNode.getChild(1));
		else
			nLineages += activeLineagesPerTransmissionTreeNode.get(trNode.getChild(1));

		return nLineages;
	}

	Node getMultifurcationParent(Node child, Node parent) {
		Node multParent;
		if (!parent.isRoot() && child.getHeight() == parent.getParent().getHeight())
			multParent = getMultifurcationParent(child, parent.getParent());
		else
			multParent = parent;

		return multParent;
	}

	@Override
	protected boolean requiresRecalculation() {
		eventListDirty = true;
		return true;
	}

	@Override
	protected void restore() {
		List<GeneTreeEvent> tmp = geneTreeEventList;
		geneTreeEventList = storedGeneTreeEventList;
		storedGeneTreeEventList = tmp;

		super.restore();
	}

	@Override
	protected void store() {
		storedGeneTreeEventList.clear();
		storedGeneTreeEventList.addAll(geneTreeEventList);

		super.store();
	}

   
}