package mcmc.gibbs

import breeze.linalg.{DenseMatrix, DenseVector}
import structure.DVStructure


case class FullState(acoefs: DenseVector[Double], bcoefs: DenseVector[Double], thcoefs: DenseMatrix[Double], indics: DenseMatrix[Double],finalCoefs: DenseMatrix[Double], mt: DenseVector[Double], tauabth: DenseVector[Double])

case class FullStateList(fstateL: List[FullState])