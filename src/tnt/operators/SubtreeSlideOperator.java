/*
 * Copyright (C) 2019. Tim Vaughan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

package tnt.operators;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import pitchfork.Pitchforks;
import starbeast2.SpeciesTreeInterface;
import tnt.distribution.GeneTreeIntervals;

@Description("Implements a version of BEAST's subtree slide operator which " +
        "is applicable to trees with hard polytomies.")
public class SubtreeSlideOperator extends TreeOperator {

    public Input<Double> relSizeInput = new Input<>("relSize",
            "Size of slide window, relative to tree height.",
            0.15);

    public Input<Double> probCoalAttachInput = new Input<>("probCoalAttach",
            "Probability of attaching to the nearest coalescent node following slide.",
            0.1);

	public Input<Double> probMultiMergerInput = new Input<>("probMultiMerger",
			"Probability of attaching to the same time as the closest(in time) coalescent node following slide.",
			0.1);

	public Input<GeneTreeIntervals> geneTreeIntervalsInput = new Input<>("geneTreeIntervals",
			"intervals for a gene tree", Validate.REQUIRED);

    Tree tree;
	double probCoalAttach, relSize, probMultiMerger;
	double[] probArray;
	GeneTreeIntervals intervals;
	HashMap<Integer, Integer> geneTreeNodeAssignment;
	List<Node> geneTreeInternalNodesList;
	SpeciesTreeInterface transmissionTree;
	Node multiMergerSibling;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        probCoalAttach = probCoalAttachInput.get();
		probMultiMerger = probMultiMergerInput.get();
		if (probCoalAttach + probMultiMerger > 1) {
			System.err.print("Sum of probCoalAttach and probMultiMerger must be less or equal to one.");
			System.exit(1);
		}
		probArray = new double[] { probCoalAttach, probMultiMerger, 1.0 - probCoalAttach - probMultiMerger };
        relSize = relSizeInput.get();

		intervals = geneTreeIntervalsInput.get();
    }

    @Override
    public double proposal() {

        double logHR = 0.0;

//        System.out.println("Count: " + (++count));

        // Select base node of edge to move:

        List<Node> logicalNodes = Pitchforks.getTrueNodes(tree);
        Node edgeBaseNode;
        do {
            edgeBaseNode = logicalNodes.get(Randomizer.nextInt(logicalNodes.size()));
        } while (edgeBaseNode.isRoot());

        // Forward HR contribution of edge node selection:

        logHR -= Math.log(1.0/(logicalNodes.size()-1.0));

        // Slide edge in randomly chosen direction:

        boolean isSlideUp = Randomizer.nextBoolean();

		geneTreeNodeAssignment = intervals.getGeneTreeNodeAssignment();
		geneTreeInternalNodesList = tree.getInternalNodes();
		geneTreeInternalNodesList = geneTreeInternalNodesList.stream()
				.sorted(Comparator.comparingDouble(n -> n.getHeight()))
				.collect(Collectors.toList());

        if (isSlideUp)
            logHR += slideUp(edgeBaseNode);
        else
            logHR += slideDown(edgeBaseNode);

        // Reverse HR contribution of edge node selection:

        logHR += Math.log(1.0/(Pitchforks.getTrueNodes(tree).size()-1.0));

        return logHR;

    }

    double getCurrentLambda() {
        return 1.0/(relSize*tree.getRoot().getHeight());
    }

    private double slideUp(Node edgeBaseNode) {

        Node edgeParentNode = edgeBaseNode.getParent();
        Node edgeSisterNode = getOtherChild(edgeParentNode, edgeBaseNode);

		AttachmentPoint newAttachmentPoint;
		try {
			newAttachmentPoint = getOlderAttachmentPoint(Pitchforks.getLogicalNode(edgeParentNode));
		} catch (AttachmentException ex) {
			return Double.NEGATIVE_INFINITY;
		}

        // Record old attachment point

        /* (Complexity is due to convention that polytomy attachments are to
            the beginning, not the end, of the attachment edge.) */

        AttachmentPoint oldAttachmentPoint = new AttachmentPoint();
        oldAttachmentPoint.attachmentHeight = edgeParentNode.getHeight();
        Node edgeLogicalParentNode = Pitchforks.getLogicalNode(edgeParentNode);
        if (Pitchforks.isPolytomy(edgeParentNode) && edgeLogicalParentNode != edgeParentNode)
            oldAttachmentPoint.attachmentEdgeBase = Pitchforks.getLogicalNode(edgeParentNode);
        else
            oldAttachmentPoint.attachmentEdgeBase = edgeSisterNode;

        // Topology modification:

        if (edgeParentNode != newAttachmentPoint.attachmentEdgeBase) {

            Node grandParent = edgeParentNode.getParent();
            grandParent.removeChild(edgeParentNode);
            edgeParentNode.removeChild(edgeSisterNode);
            grandParent.addChild(edgeSisterNode);
            edgeParentNode.setParent(null);

            if (!newAttachmentPoint.attachmentEdgeBase.isRoot()) {
                Node newGrandParent = newAttachmentPoint.attachmentEdgeBase.getParent();
                newGrandParent.removeChild(newAttachmentPoint.attachmentEdgeBase);
                newGrandParent.addChild(edgeParentNode);
            }
            edgeParentNode.addChild(newAttachmentPoint.attachmentEdgeBase);

            if (edgeParentNode.isRoot())
                tree.setRoot(edgeParentNode);
        }
        edgeParentNode.setHeight(newAttachmentPoint.attachmentHeight);

        // Probability of reverse move:

        computeYoungerAttachmentPointProb(oldAttachmentPoint, edgeParentNode);

        return oldAttachmentPoint.logProb - newAttachmentPoint.logProb;
    }

    private double slideDown(Node edgeBaseNode) {
        Node edgeParentNode = edgeBaseNode.getParent();
        Node edgeSisterNode = getOtherChild(edgeParentNode, edgeBaseNode);

        AttachmentPoint newAttachmentPoint;
        try {
            newAttachmentPoint = getYoungerAttachmentPoint(edgeBaseNode,
                    Pitchforks.getLogicalNode(edgeParentNode));
        } catch (AttachmentException ex) {
            return Double.NEGATIVE_INFINITY;
        }

        AttachmentPoint oldAttachmentPoint = new AttachmentPoint();
        oldAttachmentPoint.attachmentHeight = edgeParentNode.getHeight();
        oldAttachmentPoint.attachmentEdgeBase = Pitchforks.getLogicalNode(edgeSisterNode);

        // Topology modification

        if (edgeSisterNode != newAttachmentPoint.attachmentEdgeBase) {
            if (!edgeParentNode.isRoot()) {
                Node grandParent = edgeParentNode.getParent();
                grandParent.removeChild(edgeParentNode);
                edgeParentNode.removeChild(edgeSisterNode);
                grandParent.addChild(edgeSisterNode);
                edgeParentNode.setParent(null);
            } else {
                edgeParentNode.removeChild(edgeSisterNode);
                edgeSisterNode.setParent(null);
            }

            Node newGrandParent = newAttachmentPoint.attachmentEdgeBase.getParent();
            newGrandParent.removeChild(newAttachmentPoint.attachmentEdgeBase);
            newGrandParent.addChild(edgeParentNode);
            edgeParentNode.addChild(newAttachmentPoint.attachmentEdgeBase);

            if (edgeSisterNode.isRoot())
                tree.setRoot(edgeSisterNode);
        } else {
            // If topology is unchanged, node below edge supporting original
            // attachment will be the original edge parent node:

            oldAttachmentPoint.attachmentEdgeBase = edgeParentNode;
        }
        edgeParentNode.setHeight(newAttachmentPoint.attachmentHeight);

        // Probability of reverse move

        computeOlderAttachmentPointProb(oldAttachmentPoint, edgeParentNode);

        return oldAttachmentPoint.logProb - newAttachmentPoint.logProb;
    }


    static class AttachmentPoint {
        Node attachmentEdgeBase;
        double attachmentHeight;
        double logProb = 0;

        @Override
        public String toString() {
            return "attachmentEdgeBase: " + attachmentEdgeBase.getNr() + ", " +
                    "attachmentHeight: " + attachmentHeight + ", " +
                    "logProb: " + logProb;
        }
    }

    static class AttachmentException extends Exception { }

	AttachmentPoint getOlderAttachmentPoint(Node startNode) throws AttachmentException {

        double lambda = getCurrentLambda();

        AttachmentPoint ap = new AttachmentPoint();

        ap.attachmentEdgeBase = startNode;
        while(true) {
            Node logicalParent = Pitchforks.getLogicalParent(ap.attachmentEdgeBase);

			if (logicalParent != null) {
				int s = Randomizer.randomChoicePDF(probArray);

				if (s == 0) {
					ap.attachmentEdgeBase = logicalParent;
					ap.attachmentHeight = logicalParent.getHeight();
					ap.logProb += Math.log(probCoalAttach);
					break;
				} else if (s == 1) {
					Node trNode = intervals.transmissionTreeInput.get()
							.getNode(geneTreeNodeAssignment.get(ap.attachmentEdgeBase.getNr()));
					List<Node> possibleTrNodes = new ArrayList<Node>();
					Node trParent = trNode;
					possibleTrNodes.add(trParent);
					while (!trParent.isRoot()) {
						trParent = trParent.getParent();
						possibleTrNodes.add(trParent);
					}
					multiMergerSibling = findMultiMergerCandidate(ap.attachmentEdgeBase, possibleTrNodes, true);
					if (multiMergerSibling == null)
						throw new AttachmentException();

					Node tmpParent = Pitchforks.getLogicalParent(ap.attachmentEdgeBase);
					Node tmpBase = ap.attachmentEdgeBase;
					
					do {
						if (tmpParent.getHeight() < multiMergerSibling.getHeight()) {
							tmpBase = tmpParent;
							tmpParent = tmpBase.getParent();
						} else
							break;

					} while (!tmpBase.isRoot());
					ap.attachmentEdgeBase = tmpBase;
					ap.attachmentHeight = multiMergerSibling.getHeight();
					ap.logProb += Math.log(probMultiMerger);
					break;
				} else {
					ap.logProb += Math.log(probArray[s]);
				}
			}

//            if (logicalParent != null) {
//                if (Randomizer.nextDouble() < probCoalAttach+probMultiMerger) {
//
//                    ap.attachmentEdgeBase = logicalParent;
//                    ap.attachmentHeight = logicalParent.getHeight();
//                    ap.logProb += Math.log(probCoalAttach);
//                    break;
//                } else {
//                    ap.logProb += Math.log(1-probCoalAttach);
//                }
//            }

            double delta = Randomizer.nextExponential(lambda);

            if (logicalParent == null || delta < ap.attachmentEdgeBase.getLength()) {
                ap.logProb += -lambda*delta + Math.log(lambda);
				ap.attachmentHeight = ap.attachmentEdgeBase.getHeight() + delta; // we need to allow delta be zero?
                break;
            }

            ap.logProb += -lambda*ap.attachmentEdgeBase.getLength();
            ap.attachmentEdgeBase = logicalParent;
        }

        return ap;
    }

	private Node findMultiMergerCandidate(Node attachmentEdgeBase, List<Node> trNodeList, boolean up) {
		int idx = geneTreeInternalNodesList.indexOf(attachmentEdgeBase);
		int start = up ? idx + 1 : 0;
		int end = up ? geneTreeInternalNodesList.size() : idx;

		for (int i = start; i < end; i++) {
			Node geneNode = Pitchforks.getLogicalNode(geneTreeInternalNodesList.get(i));
			if (geneNode.isFake() || geneNode.isLeaf())
				continue;
			Node trNode = intervals.transmissionTreeInput.get()
					.getNode(geneTreeNodeAssignment.get(geneNode.getNr()));

			if (up && geneNode.getAllChildNodesAndSelf().indexOf(attachmentEdgeBase) == -1
					&& trNodeList.indexOf(trNode) != -1) {
				return geneNode;
			} else if (!up && attachmentEdgeBase.getAllChildNodesAndSelf().indexOf(geneNode) == -1
					&& trNodeList.indexOf(trNode) != -1) {
				return geneNode;
			}
		}
		return null;
	}

	void computeOlderAttachmentPointProb(AttachmentPoint ap, Node startNode) {

        double lambda = getCurrentLambda();

        ap.logProb = 0;

        Node currentEdgeBase = startNode;
        Node logicalParent;
        while(true) {
            logicalParent = Pitchforks.getLogicalParent(currentEdgeBase);

            if (logicalParent != null) {
				if (ap.attachmentHeight == multiMergerSibling.getHeight()) {
					ap.logProb += Math.log(probMultiMerger);
				} else if (ap.attachmentHeight <= logicalParent.getHeight()) {

                    if (ap.attachmentHeight == logicalParent.getHeight()) {
                        ap.logProb += Math.log(probCoalAttach);
//					} else if (ap.attachmentHeight == multiMergerSibling.getHeight()) {
//						ap.logProb += Math.log(probMultiMerger);
					} else {
						ap.logProb += Math.log(1 - probCoalAttach - probMultiMerger)
                                -lambda*(ap.attachmentHeight - currentEdgeBase.getHeight())
                                + Math.log(lambda);
                    }
                    break;

//				} else if (ap.attachmentHeight == multiMergerSibling.getHeight()) {
//					ap.logProb += Math.log(probMultiMerger);
				} else {
					ap.logProb += Math.log(1.0 - probCoalAttach - probMultiMerger)
                            -lambda*currentEdgeBase.getLength();
                }
            } else {
                ap.logProb += -lambda*(ap.attachmentHeight - currentEdgeBase.getHeight())
                        + Math.log(lambda);
                break;
            }

            currentEdgeBase = logicalParent;
        }
    }

    AttachmentPoint getYoungerAttachmentPoint(Node edgeBaseNode,
                                                      Node startNode) throws AttachmentException {

        double lambda = getCurrentLambda();

        AttachmentPoint ap = new AttachmentPoint();

        ap.attachmentEdgeBase = startNode;
        while(true) {
            List<Node> logicalChildren = Pitchforks.getLogicalChildren(ap.attachmentEdgeBase);
            if (ap.attachmentEdgeBase == startNode)
                logicalChildren.remove(edgeBaseNode);

            ap.attachmentEdgeBase = Pitchforks.randomChoice(logicalChildren);
            ap.logProb += Math.log(1.0/logicalChildren.size());

            if (!ap.attachmentEdgeBase.isLeaf()) {
				int s = Randomizer.randomChoicePDF(probArray);

				if (s == 0) {
					ap.attachmentHeight = ap.attachmentEdgeBase.getHeight();
					ap.logProb += Math.log(probArray[s]);
					break;
				}
				if (s == 1) {
					Node trNode = intervals.transmissionTreeInput.get()
							.getNode(geneTreeNodeAssignment.get(ap.attachmentEdgeBase.getNr()));
					List<Node> possibleTrNodes = elidgibleTrNodesDown(trNode);
					multiMergerSibling = findMultiMergerCandidate(ap.attachmentEdgeBase, possibleTrNodes, false);
					if (multiMergerSibling == null)
						throw new AttachmentException();

					List<Node> tmpLogicalChildren = Pitchforks.getLogicalChildren(ap.attachmentEdgeBase);
					ap.attachmentEdgeBase = Pitchforks.randomChoice(tmpLogicalChildren);
		            ap.logProb += Math.log(1.0/logicalChildren.size());

					do {
						if (ap.attachmentEdgeBase.getHeight() > multiMergerSibling.getHeight()) {
							if (ap.attachmentEdgeBase.isLeaf())
								throw new AttachmentException();
							tmpLogicalChildren = Pitchforks.getLogicalChildren(ap.attachmentEdgeBase);
							if (tmpLogicalChildren.size() == 0)
								System.out.println();
							ap.attachmentEdgeBase = Pitchforks.randomChoice(tmpLogicalChildren);
							ap.logProb += Math.log(1.0 / logicalChildren.size());

						} else
							break;

					} while (!ap.attachmentEdgeBase.isLeaf());

				} else {
					ap.logProb += Math.log(probArray[s]);
				}
//
//                if (Randomizer.nextDouble() < probCoalAttach) {
//                    ap.attachmentHeight = ap.attachmentEdgeBase.getHeight();
//                    ap.logProb += Math.log(probCoalAttach);
//                    break;
//                } else {
//                    ap.logProb += Math.log(1-probCoalAttach);
//                }
			} else if (Randomizer.nextDouble() < probMultiMerger) {
				Node trNode = intervals.transmissionTreeInput.get()
						.getNode(geneTreeNodeAssignment.get(ap.attachmentEdgeBase.getNr()));
				List<Node> possibleTrNodes = elidgibleTrNodesDown(trNode);
				multiMergerSibling = findMultiMergerCandidate(ap.attachmentEdgeBase, possibleTrNodes, false);
				if (multiMergerSibling == null)
					throw new AttachmentException();

				List<Node> tmpLogicalChildren = Pitchforks.getLogicalChildren(ap.attachmentEdgeBase);
				ap.attachmentEdgeBase = Pitchforks.randomChoice(tmpLogicalChildren);
				ap.logProb += Math.log(1.0 / logicalChildren.size());

				while (true) {
					if (ap.attachmentEdgeBase.getHeight() > multiMergerSibling.getHeight()) {
						tmpLogicalChildren = Pitchforks.getLogicalChildren(ap.attachmentEdgeBase);
						ap.attachmentEdgeBase = Pitchforks.randomChoice(tmpLogicalChildren);
						ap.logProb += Math.log(1.0 / logicalChildren.size());

					} else
						break;
				}
			} else {
				ap.logProb += Math.log(1 - probMultiMerger);
			}

            double delta = Randomizer.nextExponential(lambda);

            if (delta < ap.attachmentEdgeBase.getLength()) {
                ap.logProb += -lambda*delta + Math.log(lambda);
                ap.attachmentHeight = ap.attachmentEdgeBase.getHeight()
                        + (ap.attachmentEdgeBase.getLength() - delta);
                break;
            }

            if (ap.attachmentEdgeBase.isLeaf())
                throw new AttachmentException();

            ap.logProb += -lambda*ap.attachmentEdgeBase.getLength();
        }

        if (ap.attachmentHeight < edgeBaseNode.getHeight())
            throw new AttachmentException();

        return ap;
    }

	List<Node> elidgibleTrNodesDown(Node trNode){
    	List<Node> possibleTrNodes = new ArrayList<Node>(); 
    	possibleTrNodes.add(trNode);
		
		if (!trNode.isLeaf()) {
			Node trChild1 = trNode.getChild(0);
			possibleTrNodes.addAll(elidgibleTrNodesDown(trChild1));

			Node trChild2 = trNode.getChild(1);
			if (trChild2 != null)
				possibleTrNodes.addAll(elidgibleTrNodesDown(trChild2));
		}
    	
		return possibleTrNodes;
    }

    void computeYoungerAttachmentPointProb(AttachmentPoint ap,
                                                   Node startNode) {

        double lambda = getCurrentLambda();

        ap.logProb = 0.0;

        Node currentEdgeBase = ap.attachmentEdgeBase;

        do {
            if (currentEdgeBase == null)
                throw new IllegalStateException("Probability calculation loop failed to find startNode.");

			if (multiMergerSibling != null && ap.attachmentHeight == multiMergerSibling.getHeight()) {
				ap.logProb += Math.log(probMultiMerger);
			} else if (currentEdgeBase.getHeight() <= ap.attachmentHeight) {

//				if (ap.attachmentHeight == multiMergerSibling.getHeight()) {
//					ap.logProb += Math.log(probMultiMerger);
//				} else 
				if (ap.attachmentHeight > currentEdgeBase.getHeight()) {
                    ap.logProb += -lambda * (currentEdgeBase.getParent().getHeight() - ap.attachmentHeight)
                            + Math.log(lambda);
                    if (!currentEdgeBase.isLeaf())
						ap.logProb += Math.log(1.0 - probCoalAttach - probMultiMerger);
                } else
                    ap.logProb += Math.log(probCoalAttach);

            } else {
                ap.logProb += Math.log(1.0 - probCoalAttach)
                        - lambda*currentEdgeBase.getLength();
            }

			if (currentEdgeBase.isRoot())
				System.out.println();

            currentEdgeBase = Pitchforks.getLogicalParent(currentEdgeBase);
			if (currentEdgeBase == null || currentEdgeBase.isLeaf())
				System.out.println();

            List<Node> logicalChildren = Pitchforks.getLogicalChildren(currentEdgeBase);

            if (currentEdgeBase == startNode)
                ap.logProb += Math.log(1.0 / (logicalChildren.size() - 1));
            else
                ap.logProb += Math.log(1.0 / logicalChildren.size());

        } while (currentEdgeBase != startNode);
    }


}
