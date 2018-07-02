/*
 * OldAbstractCoalescentLikelihood.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package beast.evolution.tree.coalescent;


import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeDistribution;
import beast.math.Binomial;
import beast.util.ComparableDouble;
import beast.util.HeapSort;

import java.util.ArrayList;

import static beast.evolution.tree.coalescent.IntervalType.COALESCENT;
import static beast.evolution.tree.coalescent.IntervalType.MIGRATION;
import static beast.evolution.tree.coalescent.IntervalType.NOTHING;

/**
 * Forms a base class for a number of coalescent likelihood calculators.
 * <p/>
 * It is imported from BEAST 1 where it is used as a base
 * by a number of other classes (i.e., BayesianSkylineLikelihood, GMRFSkyrideLikelihood).
 * <p/>
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: CoalescentLikelihood.java,v 1.43 2006/07/28 11:27:32 rambaut Exp $
 */
public class OldAbstractCoalescentLikelihood extends TreeDistribution {

    // PUBLIC STUFF

//    protected MultiLociTreeSet treesSet = null;
    protected Tree tree;
    protected TreeIntervals treeIntervals;

    private boolean storedIntervalsKnown = false;

    private PopulationFunction demoFunction = null;

//    public enum CoalescentEventType {
//        /**
//         * Denotes an interval after which a coalescent event is observed
//         * (i.e. the number of lineages is smaller in the next interval)
//         */
//        COALESCENT,
//        /**
//         * Denotes an interval at the end of which a new sample addition is
//         * observed (i.e. the number of lineages is larger in the next interval).
//         */
//        NEW_SAMPLE,
//        /**
//         * Denotes an interval at the end of which nothing is
//         * observed (i.e. the number of lineages is the same in the next interval).
//         */
//        NOTHING
//    }


    public OldAbstractCoalescentLikelihood() {
        treeInput.setRule(Input.Validate.FORBIDDEN);
    }

    @Override
    public void initAndValidate() {
        if (treeInput.get() != null) {
            throw new IllegalArgumentException("only tree treeIntervals (not tree) be specified");
        }
        treeIntervals = treeIntervalsInput.get();
    }

    // **************************************************************
    // Extendable methods
    // **************************************************************

    /**
     * @param tree given tree
     * @return the node ref of the MRCA of this coalescent prior in the given tree (i.e. root of tree)
     */
    public Node getMRCAOfCoalescent(Tree tree) {
        return tree.getRoot();
    }

    /**
     * @param tree given tree
     * @return an array of noderefs that represent the MRCAs of subtrees to exclude from coalescent prior.
     *         May return null if no subtrees should be excluded.
     */
    public Node[] getExcludedMRCAs(Tree tree) {
        return null;
    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    /**
     * Stores the precalculated state: in this case the intervals
     */
    public void store() {
        if (tree != null) {
            System.arraycopy(treeIntervals.intervals, 0, treeIntervals.storedIntervals, 0, treeIntervals.intervals.length);
            System.arraycopy(treeIntervals.lineageCounts, 0, treeIntervals.storedLineageCounts, 0, treeIntervals.lineageCounts.length);
            storedIntervalsKnown = treeIntervals.intervalsKnown;
            treeIntervals.storedIntervalCount = treeIntervals.intervalCount;
//            storedLikelihoodKnown = likelihoodKnown;
        }
//        else if (treesSet != null) {
//            treesSet.storeTheState();
//        }
        super.store();
    }

    /**
     * Restores the precalculated state: that is the intervals of the tree.
     */
    public void restore() {
        if (tree != null) {
            System.arraycopy(treeIntervals.storedIntervals, 0, treeIntervals.intervals, 0, treeIntervals.storedIntervals.length);
            System.arraycopy(treeIntervals.storedLineageCounts, 0, treeIntervals.lineageCounts, 0, treeIntervals.storedLineageCounts.length);
            treeIntervals.intervalsKnown = storedIntervalsKnown;
            treeIntervals.intervalCount = treeIntervals.storedIntervalCount;
        }
//        else if (treesSet != null) {
//            treesSet.restoreTheState();
//        }

//        likelihoodKnown = storedLikelihoodKnown;

//        if (!treeIntervals.intervalsKnown) {
//            likelihoodKnown = false;
//        }
        super.restore();
    }

//    protected final void acceptState() {
//    } // nothing to do

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

//    public double getLogLikelihood() {
//        if (!likelihoodKnown) {
//            logP = calculateLogLikelihood();
//            likelihoodKnown = true;
//        }
//        return logP;
//    }

//    public final void makeDirty() {
//        likelihoodKnown = false;
//        treeIntervals.intervalsKnown = false;
//    }

    /**
     * @return the log likelihood of this set of coalescent intervals,
     *         given a demographic model
     */
    @Override
    public double calculateLogP() {

//        if (treesSet != null) {
//            final int nTrees = treesSet.nLoci();
//            final DemographicFunction demogFunction = demoModel.getDemographicFunction();
//            double logLike = 0.0;
//            for (int nt = 0; nt < nTrees; ++nt) {
//                final double popFactor = treesSet.getPopulationFactor(nt);
//                DemographicFunction df = popFactor != 1.0 ?
//                        new ScaledDemographic(demogFunction, popFactor) : demogFunction;
//
//                logLike += Coalescent.calculateLogLikelihood(treesSet.getTreeIntervals(nt), df);
//            }
//            return logLike;
//        }

        if (!treeIntervals.intervalsKnown) setupIntervals();

        if (demoFunction == null) return calculateAnalyticalLogLikelihood();

        logP = 0.0;

        double currentTime = 0.0;

        for (int j = 0; j < treeIntervals.intervalCount; j++) {

            logP += calculateIntervalLikelihood(demoFunction, treeIntervals.intervals[j], currentTime, treeIntervals.lineageCounts[j],
                    getIntervalType(j));

            // insert zero-length coalescent intervals
            final int diff = getCoalescentEvents(j) - 1;
            for (int k = 0; k < diff; k++) {
                logP += calculateIntervalLikelihood(demoFunction, 0.0, currentTime, treeIntervals.lineageCounts[j] - k - 1,
                        COALESCENT);
            }

            currentTime += treeIntervals.intervals[j];
        }

        return logP;
    }

    private double calculateAnalyticalLogLikelihood() {

        final double lambda = getLambda();
        final int n = tree.getLeafNodeCount();

        // assumes a 1/theta prior
        //logLikelihood = Math.log(1.0/Math.pow(lambda,n));

        // assumes a flat prior
        //double logL = Math.log(1.0/Math.pow(lambda,n-1));
        //final double logL = - Math.log(Math.pow(lambda,n-1));
        return -(n - 1) * Math.log(lambda);
    }


    /**
     * k - number of lineages
     * N - population size
     * kingsman coalescent: interval to next coalescent event x ~ exp(lambda), where lambda = C(k,2) / N
     * Like(x ; lambda) = lambda * exp(-lambda * x)
     * so Like(N) = (C(k,2)/N) * exp(- x * C(k,2)/N)
     * lg(Like(N)) = lg(C(k,2)) - lg(N) -C(k,2) * x/N
     * <p/>
     * When N changes over time N = N(t) we have lambda(t) = C(k,2)/N(t) and the likelihood equation is
     * Like(t) = lambda(t) * exp(- integral_0^t(lambda(x) dx) )
     * <p/>
     * lg(Like(t)) = -C(k,2) * integral_0^t(1/N(x) dx) + lg(C(k,2)/N(t))
     * <p/>
     * For a sample event, the likelihood is for no event until time t, and is just the first term of the above.
     *
     * @param demogFunction  the demographic function
     * @param width          the size of the coalescent interval
     * @param timeOfPrevCoal the time of previous coalescent event (going backwards in time)
     * @param lineageCount   the number of lineages spanning this coalescent interval
     * @param type           the type of coalescent event that this interval is terminated by
     * @return likelihood of a given interval,coalescent or otherwise
     */
    //TODO duplicate to BayesianSkyline.calculateIntervalLikelihood
    public static double calculateIntervalLikelihood(PopulationFunction demogFunction,
                                                     double width, double timeOfPrevCoal, int lineageCount,
                                                     IntervalType type) {
        final double timeOfThisCoal = width + timeOfPrevCoal;

        final double intervalArea = demogFunction.getIntegral(timeOfPrevCoal, timeOfThisCoal);
        final double kchoose2 = Binomial.choose2(lineageCount);
        double like = -kchoose2 * intervalArea;

        switch (type) {
            case COALESCENT:
                final double demographic = Math.log(demogFunction.getPopSize(timeOfThisCoal));
                like += -demographic;

                break;
            default:
                break;
        }

        return like;
    }



    /**
     * @return a factor lambda such that the likelihood can be expressed as
     *         1/theta^(n-1) * exp(-lambda/theta). This allows theta to be integrated
     *         out analytically. :-)
     */
    private double getLambda() {
        double lambda = 0.0;
        for (int i = 0; i < getIntervalCount(); i++) {
            lambda += (treeIntervals.intervals[i] * treeIntervals.lineageCounts[i]);
        }
        lambda /= 2;

        return lambda;
    }

    /**
     * Recalculates all the intervals from the tree model.
     * GL: made public, to give BayesianSkylineGibbsOperator access
     */
    public final void setupIntervals() {

        if (treeIntervals.intervals == null) {
            int maxIntervalCount = tree.getNodeCount();

            treeIntervals.intervals = new double[maxIntervalCount];
            treeIntervals.lineageCounts = new int[maxIntervalCount];
            treeIntervals.storedIntervals = new double[maxIntervalCount];
            treeIntervals.storedLineageCounts = new int[maxIntervalCount];
        }

        XTreeIntervals ti = new XTreeIntervals(treeIntervals.intervals, treeIntervals.lineageCounts);
        getTreeIntervals(getMRCAOfCoalescent(tree), getExcludedMRCAs(tree), ti);
        treeIntervals.intervalCount = ti.nIntervals;

        treeIntervals.intervalsKnown = true;
    }


    /**
     * Extract coalescent times and tip information into ArrayList times from tree.
     * Upon return times contain the time of each node in the subtree below top, and at the corrosponding index
     * of childs is the descendent count for that time.
     *
     * @param top          the node to start from
     * @param excludeBelow an optional array of nodes to exclude (corresponding subtrees) from density.
     * @param times        array to fill with times
     * @param childs       array to fill with descendents count
     */
    private static void collectAllTimes(Node top, Node[] excludeBelow,
                                        ArrayList<ComparableDouble> times, ArrayList<Integer> childs) {

        times.add(new ComparableDouble(top.getHeight()));
        childs.add(top.getChildCount());

        for (int i = 0; i < top.getChildCount(); i++) {
            Node child = top.getChild(i);
            if (excludeBelow == null) {
                collectAllTimes(child, excludeBelow, times, childs);
            } else {
                // check if this subtree is included in the coalescent density
                boolean include = true;
                for (Node anExcludeBelow : excludeBelow) {
                    if (anExcludeBelow.getNr() == child.getNr()) {
                        include = false;
                        break;
                    }
                }
                if (include)
                    collectAllTimes(child, excludeBelow, times, childs);
            }
        }
    }

    private class XTreeIntervals {

        public XTreeIntervals(double[] intervals, int[] lineageCounts) {
            this.intervals = intervals;
            this.lineagesCount = lineageCounts;
        }

        int nIntervals;
        final int[] lineagesCount;
        final double[] intervals;

    }

    private static void getTreeIntervals(Node root, Node[] exclude, XTreeIntervals ti) {
        double MULTIFURCATION_LIMIT = 1e-9;

        ArrayList<ComparableDouble> times = new ArrayList<ComparableDouble>();
        ArrayList<Integer> childs = new ArrayList<Integer>();
        collectAllTimes(root, exclude, times, childs);
        int[] indices = new int[times.size()];

        HeapSort.sort(times, indices);

        final double[] intervals = ti.intervals;
        final int[] lineageCounts = ti.lineagesCount;

        // start is the time of the first tip
        double start = times.get(indices[0]).doubleValue();
        int numLines = 0;
        int i = 0;
        int intervalCount = 0;
        while (i < times.size()) {

            int lineagesRemoved = 0;
            int lineagesAdded = 0;

            final double finish = times.get(indices[i]).doubleValue();
            double next = finish;

            while (Math.abs(next - finish) < MULTIFURCATION_LIMIT) {
                final int children = childs.get(indices[i]);
                if (children == 0) {
                    lineagesAdded += 1;
                } else {
                    lineagesRemoved += (children - 1);
                }
                i += 1;
                if (i == times.size()) break;

                next = times.get(indices[i]).doubleValue();
            }
            //System.out.println("time = " + finish + " removed = " + lineagesRemoved + " added = " + lineagesAdded);
            if (lineagesAdded > 0) {

                if (intervalCount > 0 || ((finish - start) > MULTIFURCATION_LIMIT)) {
                    intervals[intervalCount] = finish - start;
                    lineageCounts[intervalCount] = numLines;
                    intervalCount += 1;
                }

                start = finish;
            }
            // add sample event
            numLines += lineagesAdded;

            if (lineagesRemoved > 0) {

                intervals[intervalCount] = finish - start;
                lineageCounts[intervalCount] = numLines;
                intervalCount += 1;
                start = finish;
            }
            // coalescent event
            numLines -= lineagesRemoved;
        }

        ti.nIntervals = intervalCount;
    }

    /**
     * @return number of intervals
     */
    public final int getIntervalCount() {
        return treeIntervals.intervalCount;
    }

    /**
     * Gets an interval.
     *
     * @param i index of interval
     * @return interval length
     */
    public final double getInterval(int i) {
        if (i >= treeIntervals.intervalCount) throw new IllegalArgumentException();
        return treeIntervals.intervals[i];
    }

    /**
     * Returns the number of uncoalesced lineages within this interval.
     * Required for s-coalescents, where new lineages are added as
     * earlier samples are come across.
     *
     * @param i lineage index
     * @return number of uncoalesced lineages within this interval.
     */
    public final int getLineageCount(int i) {
        if (i >= treeIntervals.intervalCount) throw new IllegalArgumentException();
        return treeIntervals.lineageCounts[i];
    }

    /**
     * @param i interval index
     * @return the number coalescent events in an interval
     */
    public final int getCoalescentEvents(int i) {

        if (i >= treeIntervals.intervalCount) throw new IllegalArgumentException();
        if (i < treeIntervals.intervalCount - 1) {
            return treeIntervals.lineageCounts[i] - treeIntervals.lineageCounts[i + 1];
        } else {
            return treeIntervals.lineageCounts[i] - 1;
        }
    }

    /**
     * @param i interval index
     * @return the type of interval observed.
     */
    public final IntervalType getIntervalType(int i) {

        if (i >= treeIntervals.intervalCount) throw new IllegalArgumentException();
        int numEvents = getCoalescentEvents(i);

        if (numEvents > 0) return COALESCENT;
        else if (numEvents < 0) return MIGRATION;
        else return NOTHING;
    }

    /**
     * @return total height of the genealogy represented by these
     *         intervals.
     */
    public final double getTotalHeight() {

        double height = 0.0;
        for (int j = 0; j < treeIntervals.intervalCount; j++) {
            height += treeIntervals.intervals[j];
        }
        return height;
    }

    /**
     * @return whether this set of coalescent intervals is fully resolved
     *         (i.e. whether is has exactly one coalescent event in each
     *         subsequent interval)
     */
    public final boolean isBinaryCoalescent() {
        for (int i = 0; i < treeIntervals.intervalCount; i++) {
            if (getCoalescentEvents(i) != 1) return false;
        }

        return true;
    }

    /**
     * @return whether this set of coalescent intervals coalescent only
     *         (i.e. whether is has exactly one or more coalescent event in each
     *         subsequent interval)
     */
    public final boolean isCoalescentOnly() {
        for (int i = 0; i < treeIntervals.intervalCount; i++) {
            if (getCoalescentEvents(i) < 1) return false;
        }
        return true;
    }

    public String toString() {
        return getID(); // Double.toString(getLogLikelihood());
    }

}
