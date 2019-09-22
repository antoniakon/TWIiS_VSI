/**
  * Created by Antonia Kontaratou.
  * Variable selection for interaction terms. Assume that all the main effects are present, otherwise you end up with a model containing an interaction involving a variable for which there is no main effect.
  * Main + interaction effects are estimated as random effects to solve the identifiability problem.
  * Model: Xijk|mu,aj,bk,gjk,Ijk tau~N(mu+aj+bk+Ijk*gjk,tau^-1)
  */
import java.io.File
import breeze.linalg.{*, _}
import breeze.numerics._
import cats._
import breeze.stats.mean
import cats.data.State
import scala.collection.mutable.ListBuffer

object FPStateMonadVS {

  def variableSelection(noOfIter: Int, thin: Int, N: Int, SumObs: Double, structure: DVStructure, alphaPriorMean: Double,  betaPriorMean: Double, mu0: Double, tau0: Double, a: Double, b: Double, thetaPriorMean: Double, aPrior: Double, bPrior: Double, p: Double) = {

    val sampleNo = noOfIter / thin + 1 // Number of samples created from the MCMC
    val alphaLevels= structure.nj
    val betaLevels = structure.nk
    val njk = alphaLevels * betaLevels // Number of levels of interactions

    val curCount = Array(0.0)

    //Define case classes
    case class FullState(acoefs: DenseVector[Double], bcoefs: DenseVector[Double], thcoefs: DenseMatrix[Double], indics: DenseMatrix[Double],finalCoefs: DenseMatrix[Double], mt: DenseVector[Double], tauabth: DenseVector[Double])
    case class FullStateList(fstateL: List[FullState])

    // Update mu and tau
    // helper function for mu tau
    def nextmutau(oldfullState: FullState): FullState= {
      val prevtau = oldfullState.mt(1)
      val prevmu = oldfullState.mt(0)
      val varMu = 1.0 / (tau0 + N * prevtau) //the variance for mu
      val meanMu = (mu0 * tau0 + prevtau * (SumObs - sumAllMainInterEff(structure, oldfullState.acoefs, oldfullState.bcoefs, alphaLevels, betaLevels, oldfullState.thcoefs, oldfullState.indics))) * varMu
      val newmu = breeze.stats.distributions.Gaussian(meanMu, sqrt(varMu)).draw()
      val newtau = breeze.stats.distributions.Gamma(a + N / 2.0, 1.0 / (b + 0.5 * YminusMuAndEffects(structure, prevmu, oldfullState.acoefs, oldfullState.bcoefs, oldfullState.thcoefs, oldfullState.indics))).draw() //  !!!!TO SAMPLE FROM THE GAMMA DISTRIBUTION IN BREEZE THE β IS 1/β
      oldfullState.copy(mt=DenseVector(newmu,newtau))
    }

    // Update taus (taua, taub, tauInt)
    // helper function for taus
    def nexttaus(oldfullState: FullState):FullState= {
      var sumThetajk = 0.0

      def calcsumajr(n: Int): Double = {
        def go(n: Int, sumaj: Double): Double={
        if (n < 0) sumaj
        else {
          go(n - 1, sumaj + pow((oldfullState.acoefs(n) - alphaPriorMean), 2))
        }
      }
        go(n,0.0)
    }
      val sumaj = calcsumajr(alphaLevels-1)

      def calcsumbkr(n: Int): Double = {
        def go(n: Int, sumbk: Double): Double={
          if (n < 0) sumbk
          else go(n - 1, sumbk + pow((oldfullState.bcoefs(n) - betaPriorMean), 2))
        }
        go(n,0.0)
      }
      val sumbk = calcsumbkr(betaLevels-1)


      for (j <- 0 until alphaLevels) {
        for (k <- 0 until betaLevels) {
          sumThetajk = sumThetajk + pow((oldfullState.thcoefs(j, k) - thetaPriorMean), 2) // Sum used in sampling from Gamma distribution for the precision of theta/interacions
        }
      }

      val newtauAlpha = breeze.stats.distributions.Gamma(aPrior + alphaLevels / 2.0, 1.0 / (bPrior + 0.5 * sumaj)).draw() //sample the precision of alpha from gamma
      val newtauBeta = breeze.stats.distributions.Gamma(aPrior + betaLevels / 2.0, 1.0 / (bPrior + 0.5 * sumbk)).draw() // sample the precision of beta from gamma
      val newtauTheta = breeze.stats.distributions.Gamma(aPrior + njk / 2.0, 1.0 / (bPrior + 0.5 * sumThetajk)).draw() // sample the precision of the interactions gamma from gamma Distribition

      oldfullState.copy(tauabth = DenseVector(newtauAlpha, newtauBeta, newtauTheta))
    }

    // Update alpha coefficients
    // helper function for alpha coeffs
    def nextAlphaCoefs(oldfullState: FullState):FullState={
      val curAlphaEstim = (DenseVector.zeros[Double](alphaLevels))
      for (j <- 0 until alphaLevels) {
        val SXalphaj = structure.calcAlphaSum(j) // the sum of the observations that have alpha==j
        val Nj = structure.calcAlphaLength(j) // the number of the observations that have alpha==j
        val SumBeta = sumBetaEffGivenAlpha(structure, j, oldfullState.bcoefs) //the sum of the beta effects given alpha
        val SinterAlpha = sumInterEffGivenAlpha(structure, j, oldfullState.thcoefs, oldfullState.indics) //the sum of the gamma/interaction effects given alpha
        val varPalpha = 1.0 / (oldfullState.tauabth(0) + oldfullState.mt(1) * Nj) //the variance for alphaj
        val meanPalpha = (alphaPriorMean * oldfullState.tauabth(0) + oldfullState.mt(1) * (SXalphaj - Nj * oldfullState.mt(0) - SumBeta - SinterAlpha)) * varPalpha //the mean for alphaj
        curAlphaEstim(j) = breeze.stats.distributions.Gaussian(meanPalpha, sqrt(varPalpha)).draw()
      }
      oldfullState.copy(acoefs = curAlphaEstim)
    }

    // Update beta coefficients
    // helper function for beta coeffs
    def nextBetaCoefs(oldfullState: FullState):FullState={
      val curBetaEstim = (DenseVector.zeros[Double](betaLevels))
      for (k <- 0 until betaLevels) {
        val SXbetak = structure.calcBetaSum(k) // the sum of the observations that have beta==k
        val Nk = structure.calcBetaLength(k) // the number of the observations that have beta==k
        val SumAlpha = sumAlphaGivenBeta(structure, k, oldfullState.acoefs)//the sum of the alpha effects given beta
        val SinterBeta = sumInterEffGivenBeta(structure, k, oldfullState.thcoefs, oldfullState.indics) //the sum of the gamma/interaction effects given beta
        val varPbeta = 1.0 / (oldfullState.tauabth(1) + oldfullState.mt(1) * Nk) //the variance for betak
        val meanPbeta = (betaPriorMean * oldfullState.tauabth(1) + oldfullState.mt(1) * (SXbetak - Nk * oldfullState.mt(0) - SumAlpha - SinterBeta)) * varPbeta //the mean for betak
        curBetaEstim(k) = breeze.stats.distributions.Gaussian(meanPbeta, sqrt(varPbeta)).draw()
      }
      oldfullState.copy(bcoefs = curBetaEstim)
    }

    // Update indicators, interactions and final interaction coefficients
    //Helper function for indicators and interactions
    def nextIndicsInters(oldfullState: FullState):FullState= {

      val curIndicsEstim = (DenseMatrix.zeros[Double](alphaLevels,betaLevels))
      val curThetaEstim = (DenseMatrix.zeros[Double](alphaLevels,betaLevels))
      var count = 0.0
      for (j <- 0 until alphaLevels) {
        for (k <- 0 until betaLevels) {
          val Njk = structure.getDVList(j, k).length // the number of the observations that have alpha==j and beta==k
          val SXjk = structure.getDVList(j, k).sum // the sum of the observations that have alpha==j and beta==k

          val u = breeze.stats.distributions.Uniform(0, 1).draw()

          //log-sum-exp trick
          val logInitExp = oldfullState.mt(1) * oldfullState.thcoefs(j, k) * (SXjk - Njk * (oldfullState.mt(0) + oldfullState.acoefs(j) + oldfullState.bcoefs(k) + 0.5 * oldfullState.thcoefs(j, k)))
          val logProb0 = log(1.0 - p) //The log of the probability I=0
          val logProb1 = log(p) + logInitExp //The log of the probability I=1
          val maxProb = max(logProb0, logProb1) //Find the max of the two probabilities
          val scaledProb0 = exp(logProb0 - maxProb) //Scaled by subtracting the max value and exponentiating
          val scaledProb1 = exp(logProb1 - maxProb) //Scaled by subtracting the max value and exponentiating
          var newProb0 = scaledProb0 / (scaledProb0 + scaledProb1) //Normalised
          val newProb1 = scaledProb1 / (scaledProb0 + scaledProb1) //Normalised

          if (newProb0 < u) {
            //prob0: Probability for when the indicator = 0, so if prob0 < u => indicator = 1
            curIndicsEstim(j, k) = 1.0
            count += 1.0
            val varPInter = 1.0 / (oldfullState.tauabth(2) + oldfullState.mt(1) * Njk) //the variance for gammajk
            val meanPInter = (thetaPriorMean * oldfullState.tauabth(2) + oldfullState.mt(1) * (SXjk - Njk * (oldfullState.mt(0) + oldfullState.acoefs(j) + oldfullState.bcoefs(k)))) * varPInter
            curThetaEstim(j, k) = breeze.stats.distributions.Gaussian(meanPInter, sqrt(varPInter)).draw()
          }
          else {
            //Update indicator and current interactions if indicator = 0.0
            curIndicsEstim(j, k) = 0.0
            curThetaEstim(j, k) = breeze.stats.distributions.Gaussian(thetaPriorMean, sqrt(1 / oldfullState.tauabth(2))).draw() // sample from the prior of interactions
          }
        }
      }
      curCount(0)= count
      oldfullState.copy(thcoefs = curThetaEstim, indics = curIndicsEstim, finalCoefs = curThetaEstim*:*curIndicsEstim)
    }

    // Initialise case class objects
    val initmt = DenseVector[Double](0.0,1.0)
    val inittaus = DenseVector[Double](1.0,1.0,1.0)
    val initAlphaCoefs = DenseVector.zeros[Double](alphaLevels)
    val initBetaCoefs = DenseVector.zeros[Double](betaLevels)
    val initThetas = DenseMatrix.zeros[Double](alphaLevels,betaLevels)
    val initIndics = DenseMatrix.zeros[Double](alphaLevels,betaLevels)
    val initFinals = DenseMatrix.zeros[Double](alphaLevels,betaLevels)

//    def addFullStateToList(fstateList: FullStateList): FullStateList={
//      FullStateList(fstate::fstateList)
//    }
    // Calculate the new state
    @annotation.tailrec
    def calculateNewState( n:Int, fstate:FullState, fstateList:FullStateList): FullStateList={
      if (n==0) fstateList
      else{
        val latestmt = nextmutau(fstate)
        val latesttaus = nexttaus(latestmt)
        val latestalphas = nextAlphaCoefs(latesttaus)
        val latestbetas = nextBetaCoefs(latestalphas)
        val latestFullyUpdatedState = nextIndicsInters(latestbetas)
        println(n)
        //println(latestFullyUpdatedState)
        if((n % thin).equals(0)) {
          println("store")
          calculateNewState(n-1, latestFullyUpdatedState, FullStateList(latestFullyUpdatedState::fstateList.fstateL))
        }
        else calculateNewState(n-1, latestFullyUpdatedState, fstateList)
      }
    }
    calculateNewState(noOfIter, FullState(initAlphaCoefs, initBetaCoefs, initThetas, initIndics, initFinals, initmt, inittaus), FullStateList(List(FullState(initAlphaCoefs, initBetaCoefs, initThetas, initIndics, initFinals, initmt, inittaus))))
  }

  /**
    * Add all the beta effects for a given alpha.
    */
  def sumBetaEffGivenAlpha(structure: DVStructure, alphaIndex: Int, betaEff: DenseVector[Double]): Double = {
    val betaLevels = structure.nk

    def sumbga(n: Int): Double = {
      @annotation.tailrec
      def go(n:Int, sum: Double): Double={
        if (n<0) sum
        else go(n-1, sum + structure.getDVList(alphaIndex,n).length*betaEff(n))
      }
      go(n, 0.0)
    }
    val sum = sumbga(betaLevels-1)
    sum

  }

  /**
    * Add all the alpha effects for a given beta.
    */
  def sumAlphaGivenBeta(structure: DVStructure, betaIndex: Int, alphaEff: DenseVector[Double]): Double = {
    val alphaLevels = structure.nj

    def sumagb(n: Int): Double = {
      @annotation.tailrec
      def go(n:Int, sum: Double): Double={
        if (n<0) sum
        else go(n-1, sum + structure.getDVList(n,betaIndex).length*alphaEff(n))
      }
      go(n, 0.0)
    }
    val sum = sumagb(alphaLevels-1)
    sum
  }
  /**
    * Calculate the sum of all the alpha and all the beta effects for all the observations.
    */
  def sumAllMainInterEff(structure: DVStructure, alphaEff: DenseVector[Double], betaEff: DenseVector[Double], nj: Int, nk: Int, interEff: DenseMatrix[Double], indics: DenseMatrix[Double]): Double = {
    var sumInter = 0.0

    // For alpha effects
    def sumaEff(n: Int): Double = {
      @annotation.tailrec
      def go(n:Int, sum: Double): Double={
        if (n<0) sum
        else go(n-1, sum + sumAlphaGivenBeta(structure, n, alphaEff))
      }
      go(n, 0.0)
    }
    val sumAlpha = sumaEff(nk-1)

    // For beta effects
    def sumbEff(n: Int): Double = {
      @annotation.tailrec
      def go(n:Int, sum: Double): Double={
        if (n<0) sum
        else go(n-1, sum + sumBetaEffGivenAlpha(structure, n, betaEff))
      }
      go(n, 0.0)
    }
    val sumBeta = sumbEff(nj-1)

    // For Interaction effects
    for (i <- 0 until nj) {
      for (j <- 0 until nk) {
        sumInter += sumInterEff(structure, i, j, interEff, indics)
      }
    }
    sumAlpha + sumBeta + sumInter
  }

  /**
    * Add all the interaction effects for a given alpha and a given beta taking advantage of the DVStructure
    */
  def sumInterEff(structure: DVStructure, alphaIndex: Int, betaIndex: Int, interEff: DenseMatrix[Double], indics: DenseMatrix[Double]): Double = {
    val noOfElements = structure.getDVList(alphaIndex, betaIndex).length
    val sum = noOfElements*indics(alphaIndex, betaIndex) * interEff(alphaIndex, betaIndex)
    sum
  }

  /**
    * Add all the interaction effects for a given alpha.
    */
  def sumInterEffGivenAlpha(structure: DVStructure, alphaIndex: Int, interEff: DenseMatrix[Double], indics: DenseMatrix[Double]): Double = {
    val betaLevels = structure.nk

    def sumInterEffga(n: Int): Double = {
      @annotation.tailrec
      def go(n:Int, sum: Double): Double={
        if (n<0) sum
        else go(n-1, sum + structure.getDVList(alphaIndex,n).length * indics(alphaIndex, n) * interEff(alphaIndex, n))
      }
      go(n, 0.0)
    }
    val sum = sumInterEffga(betaLevels-1)
    sum
  }

  /**
    * Add all the interaction effects for a given beta.
    */
  def sumInterEffGivenBeta(structure: DVStructure, betaIndex: Int, interEff: DenseMatrix[Double], indics: DenseMatrix[Double]): Double = {
    val alphaLevels = structure.nj

    def sumInterEffag(n: Int): Double = {
      @annotation.tailrec
      def go(n:Int, sum: Double): Double={
        if (n<0) sum
        else go(n-1, sum + structure.getDVList(n,betaIndex).length * indics(n, betaIndex) * interEff(n, betaIndex))
      }
      go(n, 0.0)
    }
    val sum = sumInterEffag(alphaLevels-1)
    sum
  }

  /**
    * Calculate the Yi-mu-u_eff-n_eff- inter_effe. To be used in estimating tau
    */
  def YminusMuAndEffects(structure:DVStructure, mu: Double, alphaEff: DenseVector[Double], betaEff: DenseVector[Double], interEff: DenseMatrix[Double], indics: DenseMatrix[Double]): Double = {
    val alphaLevels = structure.nj
    val betaLevels = structure.nk
    var sum = 0.0
    for (ialpha <- 0 until alphaLevels) {
      for (ibeta <-0 until betaLevels){
        sum += structure.getDVList(ialpha, ibeta).map(x => scala.math.pow(x-mu-alphaEff(ialpha)-betaEff(ibeta)-interEff(ialpha,ibeta)*indics(ialpha,ibeta),2)).reduce((x, y) => x + y)
      }
    }
    sum
  }

  // Calculation of the execution time
  def time[A](f: => A) = {
    val s = System.nanoTime
    val ret = f
    println("time: " + (System.nanoTime - s) / 1e6 + "ms")
    ret
  }

  def main(args: Array[String]): Unit = {
    //stop execution until press enter
    readLine()

    // Read the data
    val data = csvread(new File("/home/antonia/ResultsFromCloud/Report/100919_15x20/Data/simulInter100919.csv"))
    val sampleSize = data.rows
    val y = data(::, 0)
    val sumObs = y.toArray.sum // Sum of the values of all the observations
    val alpha = data(::, 1).map(_.toInt).map(x => x - 1)
    val beta = data(::, 2).map(_.toInt).map(x => x - 1)
    val structure = new DVStructure(y, alpha, beta)

    // Parameters
    val noOfIters = 100
    val thin = 5
    val aPrior = 1
    val bPrior = 0.0001
    val alphaPriorMean = 0.0
    val betaPriorMean = 0.0
    val mu0 = 0.0
    val tau0 = 0.0001
    val a = 1
    val b = 0.0001
    val interPriorMean = 0.0 //common mean for all the interaction effects
    val p = 0.2

    val statesResults =
      time(
        variableSelection(
          noOfIters,
          thin,
          sampleSize,
          sumObs,
          structure,
          alphaPriorMean,
          betaPriorMean,
          mu0,
          tau0,
          a,
          b,
          interPriorMean,
          aPrior,
          bPrior,
          p)
      )

    val acoefficients= statesResults.fstateL.map(f=>f.acoefs)
    println(acoefficients)
    val acoefMat= DenseMatrix(acoefficients.map(_.toArray):_*)
    val meanValsAcoef = mean(acoefMat(::, *))
    println(meanValsAcoef)

//    //Find the average
//    println("alphas")
//    val acoefficients = statesResults.acoefs.acoefs
//    val acoefMat= DenseMatrix(acoefficients.map(_.toArray):_*)
//    val meanValsAcoef = mean(acoefMat(::, *))
//    println(meanValsAcoef)
//
//    println("betas")
//    val bcoefficients = statesResults.bcoefs.bcoefs
//    val bcoefMat= DenseMatrix(bcoefficients.map(_.toArray):_*)
//    val meanValsBcoef = mean(bcoefMat(::, *))
//    println(meanValsBcoef)
//
//    println("mu, tau")
//    val mtcoefficients = statesResults.mt.mt
//    val mtcoefMat= DenseMatrix(mtcoefficients.map(_.toArray):_*)
//    val meanValsmtcoef = mean(mtcoefMat(::, *))
//    println(meanValsmtcoef)
//
//    println("taus")
//    val tauscoefficients = statesResults.tauabth.tauabth
//    val tauscoefMat= DenseMatrix(tauscoefficients.map(_.toArray):_*)
//    val outputFile = new File("/home/antonia/Desktop/tausTry.csv")
//    breeze.linalg.csvwrite(outputFile, tauscoefMat, separator = ',')
//    val meanValstauscoef = mean(tauscoefMat(::, *))
//    println(meanValstauscoef)
//
//    println("thetas")
//    val thetascoefficients = statesResults.thcoefs.thcoefs.map(x=>x.toDenseVector)
//    val thetascoefMat= DenseMatrix(thetascoefficients.map(_.toArray):_*)
//    val meanValsthetascoef = mean(thetascoefMat(::, *))
//    println(meanValsthetascoef)
//
//    println("indicators")
//    val indicscoefficients = statesResults.indics.indics.map(x=>x.toDenseVector)
//    val indicscoefMat= DenseMatrix(indicscoefficients.map(_.toArray):_*)
//    val meanValsindicscoef = mean(indicscoefMat(::, *))
//    println(meanValsindicscoef)
//
//    println("finalCoefs")
//    val finalcoefficients = statesResults.finalCoefs.finalCoefs.map(x=>x.toDenseVector)
//    val finalcoefMat= DenseMatrix(finalcoefficients.map(_.toArray):_*)
//    val meanValsfinalcoef = mean(finalcoefMat(::, *))
//    println(meanValsfinalcoef)
 }

}

