package mcmc.gibbs

import java.io.{File, FileWriter, PrintWriter}
import breeze.linalg.{*, DenseMatrix, DenseVector, max}
import breeze.numerics.{exp, log, pow, sqrt}
import structure.DVStructure

/**
  * Saturated model with Gibbs sampler. Implementation for asymmetric main effects and asymmetric interactions.
  * Model: X_ijk | mu,a_j,b_k ,gamma_jk,tau  ~ N(mu + a_j + b_k +  gamma_jk , τ^−1 )
  * Using gamma priors for taua, taub, tauTheta and Normal for the main effects a_j and b_k and the effect size theta_jk
  * Asymmetric main effects: as and bs come from a different distribution
  * Asymmetric Interactions: gamma_jk is different from gamma_kj
  **/
class SaturatedAsymmetricBoth extends VariableSelection {

  override def variableSelection(info: InitialInfo) = {
    // Initialise case class objects
    val initmt = DenseVector[Double](0.0, 1.0)
    val inittaus = DenseVector[Double](1.0, 1.0, 1.0)
    val initAlphaCoefs = DenseVector.zeros[Double](info.alphaLevels)
    val initBetaCoefs = DenseVector.zeros[Double](info.betaLevels)
    val initZetaCoefs = DenseVector.zeros[Double](info.zetaLevels)
    val initThetas = DenseMatrix.zeros[Double](info.alphaLevels, info.betaLevels) //Thetas represent the interaction coefficients gamma for this case
    val initIndics = DenseMatrix.zeros[Double](info.alphaLevels, info.betaLevels)
    val initFinals = DenseMatrix.zeros[Double](info.alphaLevels, info.betaLevels)
    val loglik = 0.0

    val fullStateInit = FullState(initAlphaCoefs, initBetaCoefs, initZetaCoefs, initThetas, initIndics, initFinals, initmt, inittaus, loglik)
    calculateAllStates(info.noOfIter, info, fullStateInit)
  }

  /**
    * Function for updating mu and tau
    */
  override def nextmutau(oldfullState: FullState, info: InitialInfo): FullState = {
    val prevtau = oldfullState.mt(1)
    val varMu = 1.0 / (info.tau0 + info.N * prevtau) //the variance for mu
    val meanMu = (info.mu0 * info.tau0 + prevtau * (info.SumObs - sumAllMainInterEff(info.structure, oldfullState.acoefs, oldfullState.bcoefs, oldfullState.thcoefs))) * varMu
    val newmu = breeze.stats.distributions.Gaussian(meanMu, sqrt(varMu)).draw()
    //Use the just updated mu to estimate tau
    val newtau = breeze.stats.distributions.Gamma(info.a + info.N / 2.0, 1.0 / (info.b + 0.5 * YminusMuAndEffects(info.structure, newmu, oldfullState.acoefs, oldfullState.bcoefs, oldfullState.thcoefs))).draw() //  !!!!TO SAMPLE FROM THE GAMMA DISTRIBUTION IN BREEZE THE β IS 1/β
    oldfullState.copy(mt = DenseVector(newmu, newtau))
  }

  /**
    * Function for updating taus (taua, taub, tauInt)
    */
  override def nexttaus(oldfullState: FullState, info: InitialInfo): FullState = {
    val njk = info.structure.sizeOfStructure() // Number of levels of interactions

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
    oldfullState.thcoefs.foreachValue(thcoef => {
      sumThetajk += pow(thcoef - info.thetaPriorMean, 2) // Sum used in sampling from Gamma distribution for the precision of theta/interacions
    })
    sumThetajk -= (info.alphaLevels * info.betaLevels - njk) * pow(0 - info.betaPriorMean, 2) //For the missing effects (if any) added extra in the sum above

    val newtauAlpha = breeze.stats.distributions.Gamma(info.aPrior + info.alphaLevelsDist / 2.0, 1.0 / (info.bPrior + 0.5 * sumaj)).draw() //sample the precision of alpha from gamma
    val newtauBeta = breeze.stats.distributions.Gamma(info.aPrior + info.betaLevelsDist / 2.0, 1.0 / (info.bPrior + 0.5 * sumbk)).draw() // sample the precision of beta from gamma
    val newtauTheta = breeze.stats.distributions.Gamma(info.aPrior + njk / 2.0, 1.0 /(info.bPrior + 0.5 * sumThetajk)).draw() // sample the precision of the interactions gamma from gamma Distribition

    oldfullState.copy(tauabth = DenseVector(newtauAlpha, newtauBeta, newtauTheta))
  }

  override def nextCoefs(oldfullState: FullState, info: InitialInfo): FullState = {
    val latestAlphaCoefs = nextAlphaCoefs(oldfullState, info)
    nextBetaCoefs(latestAlphaCoefs, info)
  }

  /**
    * Function for updating alpha coefficients
    */
  def nextAlphaCoefs(oldfullState: FullState, info: InitialInfo): FullState = {
    val curAlphaEstim = DenseVector.zeros[Double](info.alphaLevels)

    info.structure.getAllItemsMappedByA().foreach(item => {
      val j = item._1
      val SXalphaj = info.structure.calcAlphaSum(j) // the sum of the observations that have alpha==j
      val Nj = info.structure.calcAlphaLength(j) // the number of the observations that have alpha==j
      val SumBeta = sumBetaEffGivenAlpha(info.structure, j, oldfullState.bcoefs) //the sum of the beta effects given alpha
      val SinterAlpha = sumInterEffGivenAlpha(info.structure, j, oldfullState.thcoefs) //the sum of the gamma effects given alpha
      val varPalpha = 1.0 / (oldfullState.tauabth(0) + oldfullState.mt(1) * Nj) //the variance for alphaj
      val meanPalpha = (info.alphaPriorMean * oldfullState.tauabth(0) + oldfullState.mt(1) * (SXalphaj - Nj * oldfullState.mt(0) - SumBeta - SinterAlpha)) * varPalpha //the mean for alphaj
      curAlphaEstim(j) = breeze.stats.distributions.Gaussian(meanPalpha, sqrt(varPalpha)).draw()
    })
    oldfullState.copy(acoefs = curAlphaEstim)
  }

  /**
    * Function for updating beta coefficients
    */
  def nextBetaCoefs(oldfullState: FullState, info: InitialInfo): FullState = {
    val curBetaEstim = DenseVector.zeros[Double](info.betaLevels)

    info.structure.getAllItemsMappedByB().foreach(item => {
      val k = item._1
      val SXbetak = info.structure.calcBetaSum(k) // the sum of the observations that have beta==k
      val Nk = info.structure.calcBetaLength(k) // the number of the observations that have beta==k
      val SumAlpha = sumAlphaGivenBeta(info.structure, k, oldfullState.acoefs) //the sum of the alpha effects given beta
      val SinterBeta = sumInterEffGivenBeta(info.structure, k, oldfullState.thcoefs) //the sum of the gamma/interaction effects given beta
      val varPbeta = 1.0 / (oldfullState.tauabth(1) + oldfullState.mt(1) * Nk) //the variance for betak
      val meanPbeta = (info.betaPriorMean * oldfullState.tauabth(1) + oldfullState.mt(1) * (SXbetak - Nk * oldfullState.mt(0) - SumAlpha - SinterBeta)) * varPbeta //the mean for betak
      curBetaEstim(k) = breeze.stats.distributions.Gaussian(meanPbeta, sqrt(varPbeta)).draw()
    })
    oldfullState.copy(bcoefs = curBetaEstim)
  }

  /**
    * Function for updating interaction coefficients
    */
  override def nextIndicsInters(oldfullState: FullState, info: InitialInfo): FullState = {

    val estimations = info.structure.map(item => ((item.a, item.b), {
      val Njk = item.list.length // the number of the observations that have alpha==j and beta==k
      val SXjk = item.list.sum // the sum of the observations that have alpha==j and beta==k

      val varPInter = 1.0 / (oldfullState.tauabth(2) + oldfullState.mt(1) * Njk) //the variance for gammajk
      val meanPInter = (info.thetaPriorMean * oldfullState.tauabth(2) + oldfullState.mt(1) * (SXjk - Njk * (oldfullState.mt(0) + oldfullState.acoefs(item.a) + oldfullState.bcoefs(item.b)))) * varPInter
      val thetaEstim = breeze.stats.distributions.Gaussian(meanPInter, sqrt(varPInter)).draw()
      thetaEstim
    }))

    val curThetaEstim = DenseMatrix.zeros[Double](info.alphaLevels, info.betaLevels)

    estimations.seq //make sure we are on single thread to access the collections above without concurrency problem
      .foreach { case (key, value) =>
        curThetaEstim(key._1, key._2) = value
      }

    oldfullState.copy(thcoefs = curThetaEstim)
  }

  /**
    * Add all the beta effects for a given alpha.
    */
  def sumBetaEffGivenAlpha(structure: DVStructure, alphaIndex: Int, betaEff: DenseVector[Double]): Double = {
    var sum = 0.0
    structure.getAllItemsForGivenA(alphaIndex).foreach(item => {
      sum += item.list.length * betaEff(item.b)
    })
    sum
  }

  /**
    * Add all the alpha effects for a given beta.
    */
  def sumAlphaGivenBeta(structure: DVStructure, betaIndex: Int, alphaEff: DenseVector[Double]): Double = {
    var sum = 0.0
    structure.getAllItemsForGivenB(betaIndex).foreach(item => {
      sum += item.list.length * alphaEff(item.a)
    })
    sum
  }

  /**
    * Calculate the sum of all the alpha and all the beta effects for all the observations.
    */
  def sumAllMainInterEff(structure: DVStructure, alphaEff: DenseVector[Double], betaEff: DenseVector[Double], interEff: DenseMatrix[Double]): Double = {
    structure.map(item => {
      item.list.length * (alphaEff(item.a) + betaEff(item.b) + interEff(item.a, item.b))
    }).sum
  }

  /**
    * Add all the interaction effects for a given alpha.
    */
  def sumInterEffGivenAlpha(structure: DVStructure, alphaIndex: Int, interEff: DenseMatrix[Double]): Double = {
    var sum = 0.0
    structure.getAllItemsForGivenA(alphaIndex).foreach(item => {
      sum += item.list.length * interEff(item.a, item.b)
    })
    sum
  }

  /**
    * Add all the interaction effects for a given beta.
    */
  def sumInterEffGivenBeta(structure: DVStructure, betaIndex: Int, interEff: DenseMatrix[Double]): Double = {
    var sum = 0.0
    structure.getAllItemsForGivenB(betaIndex).foreach(item => {
      sum += item.list.length * interEff(item.a, item.b)
    })
    sum
  }

  /**
    * Calculate the Yi-mu-u_eff-n_eff- inter_effe. To be used in estimating tau
    */
  def YminusMuAndEffects(structure: DVStructure, mu: Double, alphaEff: DenseVector[Double], betaEff: DenseVector[Double], interEff: DenseMatrix[Double]): Double = {

    structure.map(item => {
      val a = item.a
      val b = item.b
      item.list.map(x => scala.math.pow(x - mu - alphaEff(a) - betaEff(b) - interEff(a, b) , 2)).sum
    }).sum
  }

  override def getFilesDirectory(): String = "/home/antonia/ResultsFromCloud/Report/symmetricNov/asymmetricBoth"

  override def getInputFilePath(): String = getFilesDirectory.concat("/simulInterAsymmetricBoth.csv")

  override def getOutputRuntimeFilePath(): String = getFilesDirectory().concat("/ScalaRuntime1m15x20AsymmetricBothSaturated.txt")

  override def getOutputFilePath(): String = getFilesDirectory.concat("/Scala1mAsymBothSaturated15x20.csv")

  override def printTitlesToFile(info: InitialInfo): Unit = {
    val pw = new PrintWriter(new File(getOutputFilePath))

    val thetaTitles = (1 to info.betaLevels)
      .map { j => "-".concat(j.toString) }
      .map { entry =>
        (1 to info.alphaLevels).map { i => "theta".concat(i.toString).concat(entry) }.mkString(",")
      }.mkString(",")

    pw.append("mu ,tau, taua, taub, tauInt,")
      .append( (1 to info.alphaLevels).map { i => "alpha".concat(i.toString) }.mkString(",") )
      .append(",")
      .append( (1 to info.betaLevels).map { i => "beta".concat(i.toString) }.mkString(",") )
      .append(",")
      .append(thetaTitles)
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
        .append( fullstate.acoefs.toArray.map { alpha => alpha.toString }.mkString(",") )
        .append(",")
        .append( fullstate.bcoefs.toArray.map { beta => beta.toString }.mkString(",") )
        .append(",")
        .append( fullstate.thcoefs.toArray.map { theta => theta.toString }.mkString(",") )
        .append("\n")
    }
    pw.close()
  }

}
