/*
 * GMRFSkyrideBlockUpdateOperator.java
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

package beast.evolution.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.GMRFSkyrideLikelihood;
import beast.util.Randomizer;
import no.uib.cipr.matrix.*;

import java.text.DecimalFormat;
import java.util.logging.Logger;

/* A Metropolis-Hastings operator to update the log population sizes and precision parameter jointly under a Gaussian Markov random field prior
 *
 * @author Erik Bloomquist
 * @author Marc Suchard
 * @version $Id: GMRFSkylineBlockUpdateOperator.java,v 1.5 2007/03/20 11:26:49 msuchard Exp $
 */
@Description("A Metropolis-Hastings operator to update the log population sizes and precision parameter jointly under a Gaussian Markov random field prior")
public class GMRFSkyrideBlockUpdateOperator extends Operator {//AbstractCoercableOperator {
    public Input<GMRFSkyrideLikelihood> GMRFSkyrideLikelihoodInput = new Input<GMRFSkyrideLikelihood>("likelihood","GMRF Multilocus Skyride Likelihood", Input.Validate.REQUIRED);
    public Input<Double> scaleFactorInput = new Input<Double>("scaleFactor","scale factor: larger means more bold proposals", Input.Validate.REQUIRED);
    public Input<Integer> maxIterationsInput = new Input<Integer>("maxIterations","max number of iterations for proposal...", 200);
    public Input<Double> stopValueInput = new Input<Double>("stopValue","",0.01);
//    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);

    private double scaleFactor;
    private double lambdaScaleFactor;
    private int fieldLength;

    private int maxIterations;
    private double stopValue;

    private RealParameter popSizeParameter;
    private RealParameter precisionParameter;
    private RealParameter lambdaParameter;

    GMRFSkyrideLikelihood gmrfField;

    private double[] zeros;
//    private boolean optimize = false;

    private static boolean FAIL_SILENTLY = true;

    @Override
    public void initAndValidate() {
        gmrfField = GMRFSkyrideLikelihoodInput .get();
        popSizeParameter = gmrfField.getPopSizeParameter();
        precisionParameter = gmrfField.getPrecisionParameter();
        lambdaParameter = gmrfField.getLambdaParameter();

        this.scaleFactor = scaleFactorInput.get();
        lambdaScaleFactor = 0.0;
        fieldLength = popSizeParameter.getDimension();

        this.maxIterations = maxIterationsInput.get();
        this.stopValue = stopValueInput.get();
//        setWeight(weight);

        zeros = new double[fieldLength];

//        optimize = optimiseInput.get();
    }

    private double getNewLambda(double currentValue, double lambdaScale) {
        double a = Randomizer.nextDouble() * lambdaScale - lambdaScale / 2;
        double b = currentValue + a;
        if (b > 1)
            b = 2 - b;
        if (b < 0)
            b = -b;

        return b;
    }

    private double getNewPrecision(double currentValue, double scaleFactor) {
        double length = scaleFactor - 1 / scaleFactor;
        double returnValue;


        if (scaleFactor == 1)
            return currentValue;
        if (Randomizer.nextDouble() < length / (length + 2 * Math.log(scaleFactor))) {
            returnValue = (1 / scaleFactor + length * Randomizer.nextDouble()) * currentValue;
        } else {
            returnValue = Math.pow(scaleFactor, 2.0 * Randomizer.nextDouble() - 1) * currentValue;
        }

        return returnValue;
    }

    public DenseVector getMultiNormalMean(DenseVector CanonVector, BandCholesky Cholesky) {

        DenseVector tempValue = new DenseVector(zeros);
        DenseVector Mean = new DenseVector(zeros);

        UpperTriangBandMatrix CholeskyUpper = Cholesky.getU();

        // Assume Cholesky factorization of the precision matrix Q = LL^T

        // 1. Solve L\omega = b

        CholeskyUpper.transSolve(CanonVector, tempValue);

        // 2. Solve L^T \mu = \omega

        CholeskyUpper.solve(tempValue, Mean);

        return Mean;
    }

    public DenseVector getMultiNormal(DenseVector StandNorm, DenseVector Mean, BandCholesky Cholesky) {

        DenseVector returnValue = new DenseVector(zeros);

        UpperTriangBandMatrix CholeskyUpper = Cholesky.getU();

        // 3. Solve L^T v = z

        CholeskyUpper.solve(StandNorm, returnValue);

        // 4. Return x = \mu + v

        returnValue.add(Mean);

        return returnValue;
    }


    public static DenseVector getMultiNormal(DenseVector Mean, UpperSPDDenseMatrix Variance) {
        int length = Mean.size();
        DenseVector tempValue = new DenseVector(length);
        DenseVector returnValue = new DenseVector(length);
        UpperSPDDenseMatrix ab = Variance.copy();

        for (int i = 0; i < returnValue.size(); i++)
            tempValue.set(i, Randomizer.nextGaussian());

        DenseCholesky chol = new DenseCholesky(length, true);
        chol.factor(ab);

        UpperTriangDenseMatrix x = chol.getU();

        x.transMult(tempValue, returnValue);
        returnValue.add(Mean);
        return returnValue;
    }


    public static double logGeneralizedDeterminant(UpperTriangBandMatrix Matrix) {
        double returnValue = 0;

        for (int i = 0; i < Matrix.numColumns(); i++) {
            if (Matrix.get(i, i) > 0.0000001) {
                returnValue += Math.log(Matrix.get(i, i));
            }
        }

        return returnValue;
    }

    public DenseVector newtonRaphson(double[] data, DenseVector currentGamma, SymmTridiagMatrix proposedQ) {
        return newNewtonRaphson(data, currentGamma, proposedQ, maxIterations, stopValue);
    }

    public static DenseVector newNewtonRaphson(double[] data, DenseVector currentGamma, SymmTridiagMatrix proposedQ,
                                               int maxIterations, double stopValue) {
        DenseVector iterateGamma = currentGamma.copy();
        DenseVector tempValue = currentGamma.copy();

        int numberIterations = 0;


        while (gradient(data, iterateGamma, proposedQ).norm(Vector.Norm.Two) > stopValue) {
            try {
                jacobian(data, iterateGamma, proposedQ).solve(gradient(data, iterateGamma, proposedQ), tempValue);
            } catch (no.uib.cipr.matrix.MatrixNotSPDException e) {
                if (FAIL_SILENTLY) {
                    return null;
                }
                Logger.getLogger("dr.evomodel.coalescent.operators.GMRFSkyrideBlockUpdateOperator").fine("Newton-Raphson F");
                throw new RuntimeException("Newton-Raphson F.");
            } catch (no.uib.cipr.matrix.MatrixSingularException e) {
                if (FAIL_SILENTLY) {
                    return null;
                }
                Logger.getLogger("dr.evomodel.coalescent.operators.GMRFSkyrideBlockUpdateOperator").fine("Newton-Raphson F");
                throw new RuntimeException("Newton-Raphson F.");
            }
            iterateGamma.add(tempValue);
            numberIterations++;

            if (numberIterations > maxIterations) {
                if (FAIL_SILENTLY) {
                    return null;
                }
                Logger.getLogger("dr.evomodel.coalescent.operators.GMRFSkyrideBlockUpdateOperator").fine("Newton-Raphson F");
                throw new RuntimeException("Newton Raphson algorithm did not converge within " + maxIterations + " step to a norm less than " + stopValue + "\n" +
                        "Try starting BEAST with a more accurate initial tree.");
            }
        }

        Logger.getLogger("dr.evomodel.coalescent.operators.GMRFSkyrideBlockUpdateOperator").fine("Newton-Raphson S");
        return iterateGamma;

    }

    private static DenseVector gradient(double[] data, DenseVector value, SymmTridiagMatrix Q) {

        DenseVector returnValue = new DenseVector(value.size());
        Q.mult(value, returnValue);
        for (int i = 0; i < value.size(); i++) {
            returnValue.set(i, -returnValue.get(i) - 1 + data[i] * Math.exp(-value.get(i)));
        }
        return returnValue;
    }


    private static SPDTridiagMatrix jacobian(double[] data, DenseVector value, SymmTridiagMatrix Q) {
        SPDTridiagMatrix jacobian = new SPDTridiagMatrix(Q, true);
        for (int i = 0, n = value.size(); i < n; i++) {
            jacobian.set(i, i, jacobian.get(i, i) + Math.exp(-value.get(i)) * data[i]);
        }
        return jacobian;
    }

    public double proposal() {

        double currentPrecision = precisionParameter.getValue(0);
        double proposedPrecision = this.getNewPrecision(currentPrecision, scaleFactor);

        double currentLambda = this.lambdaParameter.getValue(0);
        double proposedLambda = this.getNewLambda(currentLambda, lambdaScaleFactor);

        precisionParameter.setValue(0, proposedPrecision);
        lambdaParameter.setValue(0, proposedLambda);

        DenseVector currentGamma = new DenseVector(gmrfField.getPopSizeParameter().getDoubleValues());
        DenseVector proposedGamma;

        SymmTridiagMatrix currentQ = gmrfField.getStoredScaledWeightMatrix(currentPrecision, currentLambda);
        SymmTridiagMatrix proposedQ = gmrfField.getScaledWeightMatrix(proposedPrecision, proposedLambda);


        double[] wNative = gmrfField.getSufficientStatistics();

        UpperSPDBandMatrix forwardQW = new UpperSPDBandMatrix(proposedQ, 1);
        UpperSPDBandMatrix backwardQW = new UpperSPDBandMatrix(currentQ, 1);

        BandCholesky forwardCholesky = new BandCholesky(wNative.length, 1, true);
        BandCholesky backwardCholesky = new BandCholesky(wNative.length, 1, true);

        DenseVector diagonal1 = new DenseVector(fieldLength);
        DenseVector diagonal2 = new DenseVector(fieldLength);
        DenseVector diagonal3 = new DenseVector(fieldLength);

        DenseVector modeForward = newtonRaphson(wNative, currentGamma, proposedQ.copy());

        if (modeForward == null) {
            // used to pass on an OperatorFailedException
            return Double.NEGATIVE_INFINITY;
        }

        for (int i = 0; i < fieldLength; i++) {
            diagonal1.set(i, wNative[i] * Math.exp(-modeForward.get(i)));
            diagonal2.set(i, modeForward.get(i) + 1);

            forwardQW.set(i, i, diagonal1.get(i) + forwardQW.get(i, i));
            diagonal1.set(i, diagonal1.get(i) * diagonal2.get(i) - 1);
        }

        forwardCholesky.factor(forwardQW.copy());

        DenseVector forwardMean = getMultiNormalMean(diagonal1, forwardCholesky);

        DenseVector stand_norm = new DenseVector(zeros);

        for (int i = 0; i < zeros.length; i++)
            stand_norm.set(i, Randomizer.nextGaussian());

        proposedGamma = getMultiNormal(stand_norm, forwardMean, forwardCholesky);


        for (int i = 0; i < fieldLength; i++)
            popSizeParameter.setValue(i, proposedGamma.get(i));

//        ((Parameter.Abstract) popSizeParameter).fireParameterChangedEvent();


        double hRatio = 0;

        diagonal1.zero();
        diagonal2.zero();
        diagonal3.zero();

        DenseVector modeBackward = newtonRaphson(wNative, proposedGamma, currentQ.copy());

        if (modeBackward == null) {
            // used to pass on an OperatorFailedException
            return Double.NEGATIVE_INFINITY;
        }

        for (int i = 0; i < fieldLength; i++) {
            diagonal1.set(i, wNative[i] * Math.exp(-modeBackward.get(i)));
            diagonal2.set(i, modeBackward.get(i) + 1);

            backwardQW.set(i, i, diagonal1.get(i) + backwardQW.get(i, i));
            diagonal1.set(i, diagonal1.get(i) * diagonal2.get(i) - 1);
        }

        backwardCholesky.factor(backwardQW.copy());

        DenseVector backwardMean = getMultiNormalMean(diagonal1, backwardCholesky);

        for (int i = 0; i < fieldLength; i++) {
            diagonal1.set(i, currentGamma.get(i) - backwardMean.get(i));
        }

        backwardQW.mult(diagonal1, diagonal3);

        // Removed 0.5 * 2
        hRatio += logGeneralizedDeterminant(backwardCholesky.getU()) - 0.5 * diagonal1.dot(diagonal3);
        hRatio -= logGeneralizedDeterminant(forwardCholesky.getU() ) - 0.5 * stand_norm.dot(stand_norm);


        return hRatio;
    }

    //MCMCOperator INTERFACE

//    public final String getOperatorName() {
//        return GMRFSkyrideBlockUpdateOperatorParser.BLOCK_UPDATE_OPERATOR;
//    }

    public double getCoercableParameter() {
//        return Math.log(scaleFactor);
        return Math.sqrt(scaleFactor - 1);
    }

    public void setCoercableParameter(double value) {
//        scaleFactor = Math.exp(value);
        scaleFactor = 1 + value * value;
    }

    public double getRawParameter() {
        return scaleFactor;
    }

    public double getScaleFactor() {
        return scaleFactor;
    }

    public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    public double getMinimumAcceptanceLevel() {
        return 0.1;
    }

    public double getMaximumAcceptanceLevel() {
        return 0.4;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return 0.20;
    }

    public double getMaximumGoodAcceptanceLevel() {
        return 0.30;
    }

    public final String getPerformanceSuggestion() {

//        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
        double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        double targetProb = getTargetAcceptanceProbability();
//        dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);
        final DecimalFormat formatter = new DecimalFormat("#.###");

        double sf = OperatorUtils.optimizeWindowSize(scaleFactor, prob, targetProb);

        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }


//  public DenseVector oldNewtonRaphson(double[] data, DenseVector currentGamma, SymmTridiagMatrix proposedQ) throws OperatorFailedException{
//  return newNewtonRaphson(data, currentGamma, proposedQ, maxIterations, stopValue);
//
//}
//
//public static DenseVector newtonRaphson(double[] data, DenseVector currentGamma, SymmTridiagMatrix proposedQ,
//int maxIterations, double stopValue) {
//
//DenseVector iterateGamma = currentGamma.copy();
//int numberIterations = 0;
//while (gradient(data, iterateGamma, proposedQ).norm(Vector.Norm.Two) > stopValue) {
//inverseJacobian(data, iterateGamma, proposedQ).multAdd(gradient(data, iterateGamma, proposedQ), iterateGamma);
//numberIterations++;
//}
//
//if (numberIterations > maxIterations)
//throw new RuntimeException("Newton Raphson algorithm did not converge within " + maxIterations + " step to a norm less than " + stopValue);
//
//return iterateGamma;
//}
//
//private static DenseMatrix inverseJacobian(double[] data, DenseVector value, SymmTridiagMatrix Q) {
//
//      SPDTridiagMatrix jacobian = new SPDTridiagMatrix(Q, true);
//      for (int i = 0; i < value.size(); i++) {
//          jacobian.set(i, i, jacobian.get(i, i) + Math.exp(-value.get(i)) * data[i]);
//      }
//
//      DenseMatrix inverseJacobian = Matrices.identity(jacobian.numRows());
//      jacobian.solve(Matrices.identity(value.size()), inverseJacobian);
//
//      return inverseJacobian;
//  }
}
