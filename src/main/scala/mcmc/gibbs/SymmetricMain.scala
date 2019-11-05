package mcmc.gibbs

import breeze.linalg.{*, DenseMatrix, DenseVector, max}
import structure.{DVStructure}
import breeze.numerics.{exp, log, pow, sqrt}
import breeze.stats.mean

class SymmetricMain extends VariableSelection {

  override def variableSelection(info: InitialInfo): FullStateList = {
    // Initialise case class objects
    val initmt = DenseVector[Double](0.0,1.0)
    val inittaus = DenseVector[Double](1.0,1.0)
    val initAlphaCoefs = DenseVector.zeros[Double](info.alphaLevels)
    val initBetaCoefs = DenseVector.zeros[Double](info.betaLevels)
    val initZetaCoefs = DenseVector.zeros[Double](info.zetaLevels)
    val initThetas = DenseMatrix.zeros[Double](info.alphaLevels, info.betaLevels)
    val initIndics = DenseMatrix.zeros[Double](info.alphaLevels, info.betaLevels)
    val initFinals = DenseMatrix.zeros[Double](info.alphaLevels,info.betaLevels)

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
      sumThetajk += pow(thcoef - info.thetaPriorMean, 2) // Sum used in sampling from Gamma distribution for the precision of theta/interacions
    })

    val njk = info.alphaLevelsDist * info.betaLevelsDist // Number of levels of interactions
    val newtauZeta = breeze.stats.distributions.Gamma(info.aPrior + info.zetaLevels / 2.0, 1.0 / (info.bPrior + 0.5 * sumzj)).draw() //sample the precision of alpha from gamma
    val newtauTheta = breeze.stats.distributions.Gamma(info.aPrior + njk / 2.0, 1.0 / (info.bPrior + 0.5 * sumThetajk)).draw() // sample the precision of the interactions gamma from gamma Distribition

    oldfullState.copy(tauabth = DenseVector(newtauZeta, newtauTheta))
  }

  override def nextCoefs(oldfullState: FullState, info: InitialInfo): FullState = {
    nextZetaCoefs(oldfullState, info)
  }
  // Update Zeta coefficients
  // helper function for Zeta coeffs
  def nextZetaCoefs(oldfullState: FullState, info: InitialInfo):FullState={
    val curZetaEstim = (DenseVector.zeros[Double](info.zetaLevels))
    (0 until info.zetaLevels).foreach( item => {
      val j = item
      val SXzetaj = info.structure.calcZetaSum(j) // the sum of the observations that have zeta on either side
      val Nj = info.structure.calcZetaLength(j) // the sum of the observations that have zeta on either side
      val SumZeta = sumRemZetaEffGivenZeta(info.structure, j, oldfullState.zcoefs) //the sum of the zeta effects of the rest for a given zeta
      val SinterZeta = sumInterEffGivenZeta(info.structure, j, oldfullState.thcoefs, oldfullState.indics) //the sum of the gamma/interaction effects given zeta
      val varPzeta = 1.0 / (oldfullState.tauabth(0) + oldfullState.mt(1) * Nj) //the variance for zetaj
      val meanPzeta = (info.alphaPriorMean * oldfullState.tauabth(0) + oldfullState.mt(1) * (SXzetaj - Nj * oldfullState.mt(0) - SumZeta - SinterZeta)) * varPzeta //the mean for zetaj
      curZetaEstim(j) = breeze.stats.distributions.Gaussian(meanPzeta, sqrt(varPzeta)).draw()
    })
    oldfullState.copy(zcoefs = curZetaEstim)
  }

  // Update indicators, interactions and final interaction coefficients
  //Helper function for indicators and interactions
  def nextIndicsInters(oldfullState: FullState, info: InitialInfo):FullState= {

    val curIndicsEstim = (DenseMatrix.zeros[Double](info.alphaLevels, info.betaLevels))
    val curThetaEstim = (DenseMatrix.zeros[Double](info.alphaLevels, info.betaLevels))
    var count = 0.0

    info.structure.foreach( item => {
      val Njk = item.list.length // the number of the observations that have alpha==j and beta==k
      val SXjk = item.list.sum // the sum of the observations that have alpha==j and beta==k

      val u = breeze.stats.distributions.Uniform(0, 1).draw()

      //log-sum-exp trick
      val thcoef = oldfullState.thcoefs(item.a, item.b)
      val logInitExp = oldfullState.mt(1) * thcoef * (SXjk - Njk * (oldfullState.mt(0) + oldfullState.zcoefs(item.a) + oldfullState.zcoefs(item.b) + 0.5 * thcoef))
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
        count += 1.0
        val varPInter = 1.0 / (oldfullState.tauabth(1) + oldfullState.mt(1) * Njk) //the variance for gammajk
        val meanPInter = (info.thetaPriorMean * oldfullState.tauabth(1) + oldfullState.mt(1) * (SXjk - Njk * (oldfullState.mt(0) + oldfullState.zcoefs(item.a) + oldfullState.zcoefs(item.b)))) * varPInter
        curThetaEstim(item.a, item.b) = breeze.stats.distributions.Gaussian(meanPInter, sqrt(varPInter)).draw()
      }
      else {
        //Update indicator and current interactions if indicator = 0.0
        curIndicsEstim(item.a,item.b) = 0.0
        curThetaEstim(item.a,item.b) = breeze.stats.distributions.Gaussian(info.thetaPriorMean, sqrt(1 / oldfullState.tauabth(1))).draw() // sample from the prior of interactions
      }
    })

    oldfullState.copy(thcoefs = curThetaEstim, indics = curIndicsEstim, finalCoefs = curThetaEstim*:*curIndicsEstim)
  }


  /**
    * For a given z it sums the effects of all the other zeta Effects
    */
  def sumRemZetaEffGivenZeta(structure: DVStructure, zetaIndex: Int, zetaEff: DenseVector[Double]): Double ={
    def notZeta(k1: Int, k2: Int): Int={
      if(k1!=zetaIndex) k1
      else k2
    }
    structure.getAllOtherZetasItemsForGivenZ(zetaIndex).map(elem => elem._2.length*zetaEff(notZeta(elem._1._1, elem._1._2))).reduce(_+_)
  }

  /**
    * Calculate the sum of all the alpha and all the beta effects for all the observations.
    */
  def sumAllMainInterEff(structure: DVStructure, zetaEff: DenseVector[Double], nz: Int, interEff: DenseMatrix[Double], indics: DenseMatrix[Double]): Double = {
    var suma = 0.0
    var sumb = 0.0
    var sumInter = 0.0

    // For zeta effects in first column
    structure.foreach(item => {
      suma = item.list.length*zetaEff(item.a)
    })

    // For zeta effects in second column
    structure.foreach(item => {
      sumb = item.list.length*zetaEff(item.b)
    })

    // Add all the interaction effects for a given alpha and a given beta taking advantage of the DVStructure
    structure.foreach( item => {
      sumInter += item.list.length * indics(item.a, item.b) * interEff(item.a, item.b)
    })

    suma + sumb + sumInter
  }


  /**
    * The sum of the interaction effects for a given zeta
    */
  def sumInterEffGivenZeta(structure: DVStructure, zetaIndex: Int, interEff: DenseMatrix[Double], indics: DenseMatrix[Double]): Double = {
    var sum = 0.0
    structure.getAllOtherZetasItemsForGivenZ(zetaIndex).map(elem => elem._2.length* indics(elem._1._1, elem._1._2) * interEff(elem._1._1, elem._1._2)).reduce(_+_)
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
    println("zetas")
    val zcoefficients = statesResults.fstateL.map(f => f.zcoefs)
    //println(acoefficients)
    val zcoefMat = DenseMatrix(zcoefficients.map(_.toArray): _*)
    val meanValsZcoef = mean(zcoefMat(::, *))
    println(meanValsZcoef)

    val matrices = calculateAndPrintCommons(statesResults)

    // Save the results to a csv file
    //    val mergedMatrix = DenseMatrix.horzcat(matrices(0), matrices(1), zcoefMat, matrices(2), matrices(3))
    //    saveToCSV(mergedMatrix, "/home/antonia/ResultsFromCloud/Report/symmetricOct/asymmetricBoth/asymmetricBothScalaRes.csv")
  }

  override def getInputFilePath(): String = "/home/antonia/ResultsFromCloud/Report/Symmetric/symmetricMain/simulInterSymmetricMain.csv"

}
