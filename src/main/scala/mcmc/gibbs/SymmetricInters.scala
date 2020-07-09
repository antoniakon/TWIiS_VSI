package mcmc.gibbs

import java.io.{File, FileWriter, PrintWriter}
import breeze.linalg.{DenseMatrix, DenseVector, max, upperTriangular}
import breeze.numerics.{exp, log, pow, sqrt}

/**
  * Variable selection with Gibbs sampler. Implementation for asymmetric main effects and symmetric interactions.
  * Extends AsymmetricBoth for the main effects and mu, tau, and implements updates for taus and interactions
  * Model: X_ijk | mu,a_j,b_k ,I_jk,theta_jk,tau  ~ N(mu + a_j + b_k + I_jk * theta_jk , τ^−1 )
  * Using gamma priors for taua, taub, tauTheta, Bermouli for the variable selection indicator I_jk and Normal for the main effects a_j and b_k and the effect size theta_jk
  * Asymmetric main effects: as and bs come from a different distribution
  * Symmetric Interactions: I_jk * theta_jk = I_kj * theta_kj
  **/
class SymmetricInters extends AsymmetricBoth {

  override def variableSelection(info: InitialInfo) = {
    // Initialise case class objects
    val initmt = DenseVector[Double](0.0, 1.0)
    val inittaus = DenseVector[Double](1.0, 1.0, 1.0)
    val initAlphaCoefs = DenseVector.zeros[Double](info.alphaLevels)
    val initBetaCoefs = DenseVector.zeros[Double](info.betaLevels)
    val initZetaCoefs = DenseVector.zeros[Double](info.zetaLevels)
    val initThetas = DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels)
    val initIndics = DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels)
    val initFinals = DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels)
    val loglik = 0.0
    val incProb = 0.0

    val fullStateInit = FullState(initAlphaCoefs, initBetaCoefs, initZetaCoefs, initThetas, initIndics, initFinals, initmt, inittaus, loglik, incProb)
    calculateAllStates(info.noOfIter, info, fullStateInit)
  }

  /**
    * Function for updating the inclusion probability p
    */
  override def nextInclusionProb(oldfullState: FullState, info: InitialInfo): FullState = {
    val njk = info.noOfInters
    val SumIjk = breeze.linalg.sum(upperTriangular(oldfullState.indics))
    val alphaP = info.pap + SumIjk
    val betaP = info.pbp + njk - SumIjk
    val newp = breeze.stats.distributions.Beta.distribution(alphaP, betaP).draw()
    oldfullState.copy(inclProb = newp)
  }

  /**
    * Function for updating taus (taua, taub, tauInt)
    */
  override def nexttaus(oldfullState: FullState, info: InitialInfo): FullState = {
    val njk = info.noOfInters // Number of levels of interactions

    var sumaj = 0.0
    oldfullState.acoefs.foreachValue(acoef => {
      sumaj += pow(acoef - info.alphaPriorMean, 2)
    })
    sumaj -= (info.alphaLevels - info.alphaLevelsDist) * pow(0 - info.alphaPriorMean, 2) //For the missing effects (if any) added extra in the sum above

    var sumbk = 0.0
    oldfullState.bcoefs.foreachValue(bcoef => {
      sumbk += pow(bcoef - info.betaPriorMean, 2)
    })
    sumbk -= (info.betaLevels - info.betaLevelsDist) * pow(0 - info.betaPriorMean, 2) //For the missing effects (if any) added extra in the sum above

    var sumThetajk = 0.0
    upperTriangular(oldfullState.thcoefs).foreachValue(thcoef => { //upperTriangular includes the main diagonal
      sumThetajk += pow(thcoef - info.thetaPriorMean, 2) // Sum used in sampling from Gamma distribution for the precision of theta/interacions
    })
    sumThetajk -= (info.noOftriangular - njk) * pow(-info.thetaPriorMean, 2) //For the missing effects (if any) added extra in the sum above

    val newtauAlpha = breeze.stats.distributions.Gamma(info.aPrior + info.alphaLevelsDist / 2.0, 1.0 / (info.bPrior + 0.5 * sumaj)).draw() //sample the precision of alpha from gamma
    val newtauBeta = breeze.stats.distributions.Gamma(info.aPrior + info.betaLevelsDist / 2.0, 1.0 / (info.bPrior + 0.5 * sumbk)).draw() // sample the precision of beta from gamma
    val newtauTheta = breeze.stats.distributions.Gamma(info.aPrior + njk / 2.0, 1.0 / (info.bPrior + 0.5 * sumThetajk)).draw() // sample the precision of the interactions gamma from gamma Distribition

    oldfullState.copy(tauabth = DenseVector(newtauAlpha, newtauBeta, newtauTheta))
  }

  /**
    * Function for updating indicators, interactions and final interaction coefficients
    */
  override def nextIndicsInters(oldfullState: FullState, info: InitialInfo): FullState = {

    val estimations = info.structureSorted.map(item => ((item.a, item.b), {
      val j = item.a
      val k = item.b

      // Number of the observations that have alpha==j and beta==k and alpha==k and beta==j
      val Njkkj = item.list.length

      // Sum of the observations that have alpha==j and beta==k and alpha==k and beta==j
      val SXjkkj = item.list.sum

      val u = breeze.stats.distributions.Uniform(0, 1).draw()

      //log-sum-exp trick
      val thcoef = oldfullState.thcoefs(item.a, item.b) // This is the same for (item.a, item.b) and (item.b, item.a) so it does not matter which one we use
      val NoOfajForbk = info.structure.calcAlphaBetaLength(j, k) //No of observations for which a==j and b==k
      val NoOfakForbj = info.structure.calcAlphaBetaLength(k, j) //No of observations for which a==k and b==j

      def returnIfExists(dv: DenseVector[Double], ind: Int) = {
        if (ind < dv.length) dv(ind)
        else 0.0
      }

      val SigmaTheta =
        if (j == k) {
          SXjkkj - Njkkj * oldfullState.mt(0) - NoOfajForbk * (returnIfExists(oldfullState.acoefs, item.a) + returnIfExists(oldfullState.bcoefs, item.b))
        } else {
          SXjkkj - Njkkj * oldfullState.mt(0) - NoOfajForbk * (returnIfExists(oldfullState.acoefs, item.a) + returnIfExists(oldfullState.bcoefs, item.b)) - NoOfakForbj * (returnIfExists(oldfullState.acoefs, item.b) + returnIfExists(oldfullState.bcoefs, item.a))
        }

      val logInitExp = oldfullState.mt(1) * thcoef * (SigmaTheta - 0.5 * Njkkj * thcoef)
      val logProb0 = log(1.0 - oldfullState.inclProb) //The log of the probability I=0
      val logProb1 = log(oldfullState.inclProb) + logInitExp //The log of the probability I=1
      val maxProb = max(logProb0, logProb1) //Find the max of the two probabilities
      val scaledProb0 = exp(logProb0 - maxProb) //Scaled by subtracting the max value and exponentiating
      val scaledProb1 = exp(logProb1 - maxProb) //Scaled by subtracting the max value and exponentiating
      var newProb0 = scaledProb0 / (scaledProb0 + scaledProb1) //Normalised
      val newProb1 = scaledProb1 / (scaledProb0 + scaledProb1) //Normalised

      if (newProb0 < u) {
        //prob0: Probability for when the indicator = 0, so if prob0 < u => indicator = 1
        val indicsEstim  = 1.0
        val varPInter = 1.0 / (oldfullState.tauabth(2) + oldfullState.mt(1) * Njkkj) //the variance for gammajk
        val meanPInter = (info.thetaPriorMean * oldfullState.tauabth(2) + oldfullState.mt(1) * SigmaTheta) * varPInter
        val thetaEstim = breeze.stats.distributions.Gaussian(meanPInter, sqrt(varPInter)).draw()
        (indicsEstim, thetaEstim)
      }
      else {
        //Update indicator and current interactions if indicator = 0.0
        val indicsEstim = 0.0
        val thetaEstim = breeze.stats.distributions.Gaussian(info.thetaPriorMean, sqrt(1 / oldfullState.tauabth(2))).draw() // sample from the prior of interactions
        (indicsEstim, thetaEstim)
      }
    }))

    val curIndicsEstim = DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels)
    val curThetaEstim = DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels)

    estimations.seq //make sure we are on single thread to access the collections above without concurrency problem
      .foreach { case (key, value) =>
        curIndicsEstim(key._1, key._2) = value._1
        curIndicsEstim(key._2, key._1) = value._1
        curThetaEstim(key._1, key._2) = value._2
        curThetaEstim(key._2, key._1) = value._2
      }

    oldfullState.copy(thcoefs = curThetaEstim, indics = curIndicsEstim, finalCoefs = curThetaEstim *:* curIndicsEstim)
  }

  override def getFilesDirectory(): String = "/home/antonia/ResultsFromCloud/Report/symmetricNov/symmetricInters"

  override def getInputFilePath(): String = getFilesDirectory.concat("/simulInterSymmetricInters.csv")

  override def getOutputRuntimeFilePath(): String = getFilesDirectory().concat("/ScalaPriorpBeta2-10SymInters1m.txt")

  override def getOutputFilePath(): String = getFilesDirectory.concat("/ScalaPriorpBeta2-10SymInters1m.csv")

  override def printTitlesToFile(info: InitialInfo): Unit = {
    val pw = new PrintWriter(new File(getOutputFilePath()))

    val thetaTitles = (1 to info.zetaLevels)
      .map { j => "-".concat(j.toString) }
      .map { entry =>
        (1 to info.zetaLevels).map { i => "theta".concat(i.toString).concat(entry) }.mkString(",")
      }.mkString(",")

    val indicsTitles = (1 to info.zetaLevels)
      .map { j => "-".concat(j.toString) }
      .map { entry =>
        (1 to info.zetaLevels).map { i => "indics".concat(i.toString).concat(entry) }.mkString(",")
      }.mkString(",")

    pw.append("mu ,tau, taua, taub, tauInt, logLik, p,")
      .append( (1 to info.alphaLevels).map { i => "alpha".concat(i.toString) }.mkString(",") )
      .append(",")
      .append( (1 to info.betaLevels).map { i => "beta".concat(i.toString) }.mkString(",") )
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
        .append(fullstate.logLik.toString)
        .append(",")
        .append(fullstate.inclProb.toString)
        .append(",")
        .append( fullstate.acoefs.toArray.map { alpha => alpha.toString }.mkString(",") )
        .append(",")
        .append( fullstate.bcoefs.toArray.map { beta => beta.toString }.mkString(",") )
        .append(",")
        .append( fullstate.finalCoefs.toArray.map { theta => theta.toString }.mkString(",") )
        .append(",")
        .append( fullstate.indics.toArray.map { ind => ind.toString }.mkString(",") )
        .append("\n")
    }
    pw.close()
  }

}
