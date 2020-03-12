package mcmc.gibbs

import java.io.{File, FileWriter, PrintWriter}
import breeze.linalg.{DenseMatrix, DenseVector, max, upperTriangular}
import breeze.numerics.{exp, log, pow, sqrt}

/**
  * Variable selection with Gibbs sampler. Implementation for symmetric main effects and symmetric interactions.
  * Extends SymmetricMain for updating z coefficients, mu and tau. Implements update for taus and interactions based on SymmetricInteractions class.
  * Model: X_ijk | mu,a_j,b_k ,I_jk,theta_jk,tau  ~ N(mu + z_j + z_k + I_jk * theta_jk , τ^−1 )
  * Using gamma priors for taua, taub, tauTheta, Bermouli for the variable selection indicator I_jk and Normal for the main effects z and the effect size theta_jk
  * Symmetric main effects: zs come from the same distribution
  * Symmetric Interactions: I_jk * theta_jk = I_kj * theta_kj
  **/
class SymmetricBoth extends SymmetricMain {

  override def variableSelection(info: InitialInfo) = {
    // Initialise case class objects
    val initmt = DenseVector[Double](0.0, 1.0)
    val inittaus = DenseVector[Double](1.0, 1.0)
    val initAlphaCoefs = DenseVector.zeros[Double](info.alphaLevels)
    val initBetaCoefs = DenseVector.zeros[Double](info.betaLevels)
    val initZetaCoefs = DenseVector.zeros[Double](info.zetaLevels)
    val initThetas = DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels)
    val initIndics = DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels)
    val initFinals = DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels)

    val fullStateInit = FullState(initAlphaCoefs, initBetaCoefs, initZetaCoefs, initThetas, initIndics, initFinals, initmt, inittaus)
    calculateAllStates(info.noOfIter, info, fullStateInit)
  }

  /**
    * Function for updating taus (tauz, tauInt)
    */
  override def nexttaus(oldfullState: FullState, info: InitialInfo): FullState = {

    //todo: check if acoef non set values create an issue
    var sumzj = 0.0
    oldfullState.zcoefs.foreachValue(zcoef => {
      sumzj += pow(zcoef - info.alphaPriorMean, 2)
    })

    //todo: check if thcoef non set values create an issue
    var sumThetajk = 0.0
    upperTriangular(oldfullState.thcoefs).foreachValue(thcoef => {
      sumThetajk += pow(thcoef - info.thetaPriorMean, 2) // Sum used in sampling from Gamma distribution for the precision of theta/interacions
    })
    sumThetajk += (info.noOfInters - info.sizeOfDouble) * pow(-info.thetaPriorMean, 2)
    val njk = info.noOfInters // Number of levels of interactions
    val newtauZeta = breeze.stats.distributions.Gamma(info.aPrior + info.zetaLevels / 2.0, 1.0 / (info.bPrior + 0.5 * sumzj)).draw() //sample the precision of alpha from gamma
    val newtauTheta = breeze.stats.distributions.Gamma(info.aPrior + njk / 2.0, 1.0 / (info.bPrior + 0.5 * sumThetajk)).draw() // sample the precision of the interactions gamma from gamma Distribition

    oldfullState.copy(tauabth = DenseVector(newtauZeta, newtauTheta))
  }

  /**
    * Function for updating indicators, interactions and final interaction coefficients
    */
  override def nextIndicsInters(oldfullState: FullState, info: InitialInfo): FullState = {

    val curIndicsEstim = (DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels))
    val curThetaEstim = (DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels))
    var count = 0.0

    info.structureSorted.foreach(item => {
      val j = item.a
      val k = item.b

      // toDo: Njkkj is the same as item.list.length from structureSorted so see maybe it is not necessary to have calcAlphaBetaLength
      // Number of the observations that have alpha==j and beta==k and alpha==k and beta==j
      val Njkkj =
      if (j == k) {
        info.structure.calcAlphaBetaLength(j, k)
      } else {
        info.structure.calcAlphaBetaLength(j, k) + info.structure.calcAlphaBetaLength(k, j)
      }

      // Sum of the observations that have alpha==j and beta==k and alpha==k and beta==j
      val SXjkkj =
        if (j == k) {
          info.structure.calcAlphaBetaSum(j, k)
        } else {
          info.structure.calcAlphaBetaSum(j, k) + info.structure.calcAlphaBetaSum(k, j)
        }

      val u = breeze.stats.distributions.Uniform(0, 1).draw()

      // log-sum-exp trick
      val thcoef = oldfullState.thcoefs(item.a, item.b) // This is the same for (item.a, item.b) and (item.b, item.a) so it does not matter which one we use
      val NoOfajForbk = info.structure.calcAlphaBetaLength(j, k) //No of observations for which a==j and b==k
      val NoOfakForbj = info.structure.calcAlphaBetaLength(k, j) //No of observations for which a==k and b==j

      def returnIfExists(dv: DenseVector[Double], ind: Int) = {
        if (ind < dv.length) dv(ind)
        else 0.0
      }

      val SigmaTheta =
        if (j == k) {
          SXjkkj - Njkkj * oldfullState.mt(0) - NoOfajForbk * (returnIfExists(oldfullState.zcoefs, item.a) + returnIfExists(oldfullState.zcoefs, item.b))
        } else {
          SXjkkj - Njkkj * oldfullState.mt(0) - NoOfajForbk * (returnIfExists(oldfullState.zcoefs, item.a) + returnIfExists(oldfullState.zcoefs, item.b)) - NoOfakForbj * (returnIfExists(oldfullState.zcoefs, item.b) + returnIfExists(oldfullState.zcoefs, item.a))
        }

      val logInitExp = oldfullState.mt(1) * thcoef * (SigmaTheta - 0.5 * Njkkj * thcoef)
      val logProb0 = log(1.0 - info.p) //The log of the probability I=0
      val logProb1 = log(info.p) + logInitExp //The log of the probability I=1
      val maxProb = max(logProb0, logProb1) //Find the max of the two probabilities
      val scaledProb0 = exp(logProb0 - maxProb) //Scaled by subtracting the max value and exponentiating
      val scaledProb1 = exp(logProb1 - maxProb) //Scaled by subtracting the max value and exponentiating
      var newProb0 = scaledProb0 / (scaledProb0 + scaledProb1) //Normalised
      val newProb1 = scaledProb1 / (scaledProb0 + scaledProb1) //Normalised

      if (newProb0 < u) {
        //prob0: Probability for when the indicator = 0, so if prob0 < u => indicator = 1
        curIndicsEstim(item.a, item.b) = 1.0
        curIndicsEstim(item.b, item.a) = curIndicsEstim(item.a, item.b)
        count += 1.0
        val varPInter = 1.0 / (oldfullState.tauabth(1) + oldfullState.mt(1) * Njkkj) //the variance for gammajk
        val meanPInter = (info.thetaPriorMean * oldfullState.tauabth(1) + oldfullState.mt(1) * SigmaTheta) * varPInter
        curThetaEstim(item.a, item.b) = breeze.stats.distributions.Gaussian(meanPInter, sqrt(varPInter)).draw()
        curThetaEstim(item.b, item.a) = curThetaEstim(item.a, item.b)
      }
      else {
        //Update indicator and current interactions if indicator = 0.0
        curIndicsEstim(item.a, item.b) = 0.0
        curIndicsEstim(item.b, item.a) = curIndicsEstim(item.a, item.b)
        curThetaEstim(item.a, item.b) = breeze.stats.distributions.Gaussian(info.thetaPriorMean, sqrt(1 / oldfullState.tauabth(1))).draw() // sample from the prior of interactions
        curThetaEstim(item.b, item.a) = curThetaEstim(item.a, item.b)
      }
    })
    oldfullState.copy(thcoefs = curThetaEstim, indics = curIndicsEstim, finalCoefs = curThetaEstim *:* curIndicsEstim)
  }

  override def getFilesDirectory(): String = "/home/antonia/ResultsFromCloud/Report/symmetricMarch/symmetricBoth"

  override def getInputFilePath(): String = getFilesDirectory.concat("/simulInterSymmetricBoth.csv")

  override def getOutputRuntimeFilePath(): String = getFilesDirectory().concat("/ScalaRuntime100kSymmetricBoth.txt")

  override def getOutputFilePath(): String = getFilesDirectory.concat("/symmetricBothScalaRes.csv")

  override def printTitlesToFile(info: InitialInfo): Unit = {
    val pw = new PrintWriter(new File(getOutputFilePath()))

    val thetaTitles = (1 to info.betaLevels)
      .map { j => "-".concat(j.toString) }
      .map { entry =>
        (1 to info.alphaLevels).map { i => "theta".concat(i.toString).concat(entry) }.mkString(",")
      }.mkString(",")

    val indicsTitles = (1 to info.betaLevels)
      .map { j => "-".concat(j.toString) }
      .map { entry =>
        (1 to info.alphaLevels).map { i => "indics".concat(i.toString).concat(entry) }.mkString(",")
      }.mkString(",")

    pw.append("mu ,tau, tauz, tauInt,")
      .append( (1 to info.zetaLevels).map { i => "zeta".concat(i.toString) }.mkString(",") )
      .append(",")
      .append(thetaTitles)
      .append(",")
      .append(indicsTitles)
      .append("\n")

    pw.close()
  }

  override def printToFile(fullStateList: FullStateList): Unit = {
    val pw = new PrintWriter(new FileWriter(getOutputFilePath(), true))

    fullStateList.fstateL.foreach { fullstate =>
      pw
        .append(fullstate.mt(0).toString)
        .append(",")
        .append(fullstate.mt(1).toString)
        .append(",")
        .append( fullstate.tauabth.toArray.map { tau => tau.toString }.mkString(",") )
        .append(",")
        .append( fullstate.zcoefs.toArray.map { alpha => alpha.toString }.mkString(",") )
        .append(",")
        .append( fullstate.thcoefs.toArray.map { theta => theta.toString }.mkString(",") )
        .append(",")
        .append( fullstate.indics.toArray.map { ind => ind.toString }.mkString(",") )
        .append("\n")
    }
    pw.close()
  }

}
