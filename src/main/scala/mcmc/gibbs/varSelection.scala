package mcmc.gibbs

import breeze.linalg.{DenseMatrix, DenseVector}
import structure.DVStructure


case class FullState(acoefs: DenseVector[Double], bcoefs: DenseVector[Double], thcoefs: DenseMatrix[Double], indics: DenseMatrix[Double],finalCoefs: DenseMatrix[Double], mt: DenseVector[Double], tauabth: DenseVector[Double])

case class FullStateList(fstateL: List[FullState])

case class InitialInfo(noOfIter: Int, thin: Int, N: Int, SumObs: Double, structure: DVStructure, alphaLevels: Int, betaLevels: Int,
                       alphaPriorMean: Double, betaPriorMean: Double, thetaPriorMean: Double, mu0: Double, tau0: Double,
                       a: Double, b: Double, aPrior: Double, bPrior: Double, p: Double)