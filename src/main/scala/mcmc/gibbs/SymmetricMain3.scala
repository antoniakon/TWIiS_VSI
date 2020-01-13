package mcmc.gibbs

import breeze.linalg.{*, DenseMatrix, DenseVector, max}
import breeze.numerics.{exp, log, pow, sqrt}
import breeze.stats.mean
import structure.DVStructure

class SymmetricMain3 extends VariableSelection {

  override def variableSelection(info: InitialInfo): FullStateList = {

    // Initialise case class objects
    val initmt = DenseVector[Double](0.0,1.0)
    val inittaus = DenseVector[Double](1.0,1.0,1.0)
    val initAlphaCoefs = DenseVector.zeros[Double](info.alphaLevels) //Not used in SymmetricMain implementation
    val initBetaCoefs = DenseVector.zeros[Double](info.betaLevels) //Not used in SymmetricMain implementation
    val initZetaCoefs = DenseVector.zeros[Double](info.zetaLevels)
    val initThetas = DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels) //zetaLevels in SymmetricMain implementation
    val initIndics = DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels) //zetaLevels in SymmetricMain implementation
    val initFinals = DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels) //zetaLevels in SymmetricMain implementation

    calculateNewState(info.noOfIter, info, FullState(initAlphaCoefs, initBetaCoefs, initZetaCoefs, initThetas, initIndics, initFinals, initmt, inittaus), FullStateList(List(FullState(initAlphaCoefs, initBetaCoefs, initZetaCoefs, initThetas, initIndics, initFinals, initmt, inittaus))))
  }

  // Update mu and tau
  // helper function for mu tau
  override def nextmutau(oldfullState: FullState, info: InitialInfo): FullState= {
    val prevtau = oldfullState.mt(1)
    val prevmu = oldfullState.mt(0)
    val varMu = 1.0 / (info.tau0 + info.N * prevtau) //the variance for mu
    val meanMu = (info.mu0 * info.tau0 + prevtau * (info.SumObs - sumAllMainInterEff(info.structure, oldfullState.zcoefs, info.zetaLevels, oldfullState.thcoefs, oldfullState.indics))) * varMu
    val newmu = breeze.stats.distributions.Gaussian(meanMu, sqrt(varMu)).draw()
    val newtau = breeze.stats.distributions.Gamma(info.a + info.N / 2.0, 1.0 / (info.b + 0.5 * YminusMuAndEffects(info.structure, prevmu, oldfullState.zcoefs, oldfullState.thcoefs, oldfullState.indics))).draw() //  !!!!TO SAMPLE FROM THE GAMMA DISTRIBUTION IN BREEZE THE β IS 1/β
    oldfullState.copy(mt=DenseVector(newmu,newtau))
  }

  // Update taus (taua, taub, tauInt)
  // helper function for taus
  override def nexttaus(oldfullState: FullState, info: InitialInfo):FullState= {

    //todo: check if acoef non set values create an issue
    var sumzj = 0.0

    oldfullState.zcoefs.foreachValue( zcoef => {
      sumzj += pow(zcoef - info.alphaPriorMean, 2)
    })

    //todo: check if thcoef non set values create an issue
    var sumThetajk = 0.0
    oldfullState.thcoefs.foreachValue(thcoef => {
      sumThetajk += pow(thcoef -info.thetaPriorMean, 2) // Sum used in sampling from Gamma distribution for the precision of theta/interacions
    })

    val njk = info.alphaLevelsDist * info.betaLevelsDist // Number of levels of interactions
    val newtauZeta = breeze.stats.distributions.Gamma(info.aPrior + info.alphaLevels / 2.0, 1.0 / (info.bPrior + 0.5 * sumzj)).draw() //sample the precision of alpha from gamma
    val newtauTheta = breeze.stats.distributions.Gamma(info.aPrior + njk / 2.0, 1.0 / (info.bPrior + 0.5 * sumThetajk)).draw() // sample the precision of the interactions gamma from gamma Distribition

    oldfullState.copy(tauabth = DenseVector(newtauZeta, newtauTheta))
  }

  override def nextCoefs(oldfullState: FullState, info: InitialInfo): FullState = {
    val latestAlphaCoefs = nextAlphaCoefs(oldfullState, info)
    nextBetaCoefs(latestAlphaCoefs, info)
  }

  // Update alpha coefficients
  // helper function for alpha coeffs
  def nextAlphaCoefs(oldfullState: FullState, info: InitialInfo):FullState={

    val curAlphaEstim = DenseVector.zeros[Double](info.alphaLevels)
    info.structure.getAllItemsMappedByA().foreach( item => {
      val j = item._1
      val SXalphaj = info.structure.calcAlphaSum(j) // the sum of the observations that have alpha==j
      val Nj = info.structure.calcAlphaLength(j) // the number of the observations that have alpha==j
      val SumBeta = sumBetaEffGivenAlpha(info.structure, j, oldfullState.bcoefs) //the sum of the beta effects given alpha
      val SinterAlpha = sumInterEffGivenAlpha(info.structure, j, oldfullState.thcoefs, oldfullState.indics) //the sum of the gamma/interaction effects given alpha
      val varPalpha = 1.0 / (oldfullState.tauabth(0) + oldfullState.mt(1) * Nj) //the variance for alphaj
      val meanPalpha = (info.alphaPriorMean * oldfullState.tauabth(0) + oldfullState.mt(1) * (SXalphaj - Nj * oldfullState.mt(0) - SumBeta - SinterAlpha)) * varPalpha //the mean for alphaj
      curAlphaEstim(j) = breeze.stats.distributions.Gaussian(meanPalpha, sqrt(varPalpha)).draw()
    })

    oldfullState.copy(acoefs = curAlphaEstim)
  }

  // Update beta coefficients
  // helper function for beta coeffs
  def nextBetaCoefs(oldfullState: FullState, info: InitialInfo):FullState={
    val curBetaEstim = DenseVector.zeros[Double](info.betaLevels)
    info.structure.getAllItemsMappedByB().foreach( item => {
      val k = item._1
      val SXbetak = info.structure.calcBetaSum(k) // the sum of the observations that have beta==k
      val Nk = info.structure.calcBetaLength(k) // the number of the observations that have beta==k
      val SumAlpha = sumAlphaGivenBeta(info.structure, k, oldfullState.acoefs)//the sum of the alpha effects given beta
      val SinterBeta = sumInterEffGivenBeta(info.structure, k, oldfullState.thcoefs, oldfullState.indics) //the sum of the gamma/interaction effects given beta
      val varPbeta = 1.0 / (oldfullState.tauabth(1) + oldfullState.mt(1) * Nk) //the variance for betak
      val meanPbeta = (info.betaPriorMean * oldfullState.tauabth(1) + oldfullState.mt(1) * (SXbetak - Nk * oldfullState.mt(0) - SumAlpha - SinterBeta)) * varPbeta //the mean for betak
      curBetaEstim(k) = breeze.stats.distributions.Gaussian(meanPbeta, sqrt(varPbeta)).draw()
    })
    oldfullState.copy(bcoefs = curBetaEstim)
  }

  // Update indicators, interactions and final interaction coefficients
  //Helper function for indicators and interactions

  override def nextIndicsInters(oldfullState: FullState, info: InitialInfo):FullState= {

    val curIndicsEstim = DenseMatrix.zeros[Double](info.alphaLevels, info.betaLevels)
    val curThetaEstim = DenseMatrix.zeros[Double](info.alphaLevels, info.betaLevels)
    var count = 0.0

    info.structure.foreach( item => {
      val Njk = item.list.length // the number of the observations that have alpha==j and beta==k
      val SXjk = item.list.sum // the sum of the observations that have alpha==j and beta==k

      val u = breeze.stats.distributions.Uniform(0, 1).draw()

      //log-sum-exp trick
      val thcoef = oldfullState.thcoefs(item.a, item.b)
      val logInitExp = oldfullState.mt(1) * thcoef * (SXjk - Njk * (oldfullState.mt(0) + oldfullState.acoefs(item.a) + oldfullState.bcoefs(item.b) + 0.5 * thcoef))
      val logProb0 = log(1.0 - info.p) //The log of the probability I=0
      val logProb1 = log(info.p) + logInitExp //The log of the probability I=1
      val maxProb = max(logProb0, logProb1) //Find the max of the two probabilities
      val scaledProb0 = exp(logProb0 - maxProb) //Scaled by subtracting the max value and exponentiating
      val scaledProb1 = exp(logProb1 - maxProb) //Scaled by subtracting the max value and exponentiating
      val newProb0 = scaledProb0 / (scaledProb0 + scaledProb1) //Normalised
//      val newProb1 = scaledProb1 / (scaledProb0 + scaledProb1) //Normalised

      if (newProb0 < u) {
        //prob0: Probability for when the indicator = 0, so if prob0 < u => indicator = 1
        curIndicsEstim(item.a, item.b) = 1.0
        count += 1.0
        val varPInter = 1.0 / (oldfullState.tauabth(2) + oldfullState.mt(1) * Njk) //the variance for gammajk
        val meanPInter = (info.thetaPriorMean * oldfullState.tauabth(2) + oldfullState.mt(1) * (SXjk - Njk * (oldfullState.mt(0) + oldfullState.acoefs(item.a) + oldfullState.bcoefs(item.b)))) * varPInter
        curThetaEstim(item.a, item.b) = breeze.stats.distributions.Gaussian(meanPInter, sqrt(varPInter)).draw()
      }
      else {
        //Update indicator and current interactions if indicator = 0.0
        curIndicsEstim(item.a,item.b) = 0.0
        curThetaEstim(item.a,item.b) = breeze.stats.distributions.Gaussian(info.thetaPriorMean, sqrt(1 / oldfullState.tauabth(2))).draw() // sample from the prior of interactions
      }
    })

    oldfullState.copy(thcoefs = curThetaEstim, indics = curIndicsEstim, finalCoefs = curThetaEstim*:*curIndicsEstim)
  }

  /**
    * Add all the beta effects for a given alpha.
    */
  def sumBetaEffGivenAlpha(structure: DVStructure, alphaIndex: Int, betaEff: DenseVector[Double]): Double = {
    var sum = 0.0
    structure.getAllItemsForGivenA(alphaIndex).foreach( item => {
      sum += item.list.length * betaEff(item.b)
    })

    sum
  }

  /**
    * Add all the alpha effects for a given beta.
    */
  def sumAlphaGivenBeta(structure: DVStructure, betaIndex: Int, alphaEff: DenseVector[Double]): Double = {
    var sum = 0.0
    structure.getAllItemsForGivenB(betaIndex).foreach( item => {
      sum += item.list.length * alphaEff(item.a)
    })

    sum
  }
  /**
    * Calculate the sum of all the zeta 1 and all the zeta 2 effects for all the observations.
    */
  def sumAllMainInterEff(structure: DVStructure, zetaEff: DenseVector[Double], nz: Int, interEff: DenseMatrix[Double], indics: DenseMatrix[Double]): Double = {
    var suma = 0.0
    var sumb = 0.0
    var sumInter = 0.0

    // For zeta effects in first column
    structure.foreach(item => {
      suma += item.list.length*zetaEff(item.a)
    })

    // For zeta effects in second column
    structure.foreach(item => {
      sumb += item.list.length*zetaEff(item.b)
    })

    // Add all the interaction effects for a given alpha and a given beta taking advantage of the DVStructure
    structure.foreach( item => {
      sumInter += item.list.length * indics(item.a, item.b) * interEff(item.a, item.b)
    })

    suma + sumb + sumInter
  }

  /**
    * Add all the interaction effects for a given alpha.
    */
  def sumInterEffGivenAlpha(structure: DVStructure, alphaIndex: Int, interEff: DenseMatrix[Double], indics: DenseMatrix[Double]): Double = {
    var sum = 0.0
    structure.getAllItemsForGivenA(alphaIndex).foreach( item => {
      sum += item.list.length * indics(item.a, item.b) * interEff(item.a, item.b)
    })

    sum
  }

  /**
    * Add all the interaction effects for a given beta.
    */
  def sumInterEffGivenBeta(structure: DVStructure, betaIndex: Int, interEff: DenseMatrix[Double], indics: DenseMatrix[Double]): Double = {
    var sum = 0.0
    structure.getAllItemsForGivenB(betaIndex).foreach( item => {
      sum += item.list.length * indics(item.a, item.b) * interEff(item.a, item.b)
    })

    sum
  }

  /**
    * Calculate the Yi-mu-u_eff-n_eff- inter_effe. To be used in estimating tau
    */
  def YminusMuAndEffects(structure:DVStructure, mu: Double, zetaEff: DenseVector[Double], interEff: DenseMatrix[Double], indics: DenseMatrix[Double]): Double = {
    var sum = 0.0

    structure.foreach( item => {
      val a = item.a
      val b = item.b
      sum += item.list.map(x => scala.math.pow(x - mu - zetaEff(a) - zetaEff(b) - interEff(a, b) * indics(a, b), 2)).sum
    })
    sum
  }
  
  override def printResults(statesResults: FullStateList): Unit = {
    println("alphas")
    val acoefficients = statesResults.fstateL.map(f => f.acoefs)
    //println(acoefficients)
    val acoefMat = DenseMatrix(acoefficients.map(_.toArray): _*)
    val meanValsAcoef = mean(acoefMat(::, *))
    println(meanValsAcoef)

    println("betas")
    val bcoefficients = statesResults.fstateL.map(f => f.bcoefs)
    val bcoefMat = DenseMatrix(bcoefficients.map(_.toArray): _*)
    val meanValsBcoef = mean(bcoefMat(::, *))
    println(meanValsBcoef)

    val matrices = calculateAndPrintCommons(statesResults)

    // Save the results to a csv file
    val mergedMatrix = DenseMatrix.horzcat(matrices(0), matrices(1), acoefMat, bcoefMat, matrices(2), matrices(3))
    saveToCSV(mergedMatrix, getFileNameToSaveResults("allCoefs"))
//    saveToCSV(matrices(0), getFileNameToSaveResults("mutau"))
//    saveToCSV(matrices(1), getFileNameToSaveResults("taus"))
//    saveToCSV(acoefMat, getFileNameToSaveResults("alphas"))
//    saveToCSV(bcoefMat, getFileNameToSaveResults("betas"))
//    saveToCSV(matrices(2), getFileNameToSaveResults("thetas"))
//    saveToCSV(matrices(3), getFileNameToSaveResults("indics"))
  }

  override protected def getFileNameToSaveResults(param: String): String = {
    val filePath = getMainFilePath.concat("/asymmetricBothCode10mScalaRes-")
    val pathToFiles = Map("mutau" -> filePath.concat("mutau.csv") ,
      "taus" -> filePath.concat("taus.csv"),
      "alphas" -> filePath.concat("alphas.csv"),
      "betas" -> filePath.concat("betas.csv"),
      "thetas" -> filePath.concat("thetas.csv"),
      "indics" -> filePath.concat("indics.csv"),
      "allCoefs" -> filePath.concat("allcoefs.csv")
    )
    pathToFiles(param)
  }

  override def getInputFilePath(): String = getMainFilePath.concat("/simulInterSymmetricMain.csv")

}
