import java.io.File

import breeze.linalg.{*, _}
import breeze.numerics._
import breeze.stats.mean

import scala.collection.mutable.ListBuffer

/**
  * Created by Antonia Kontaratou.
  * Variable selection for interaction terms. Assume that all the main effects are present, otherwise you end up with a model containing an interaction involving a variable for which there is no main effect.
  * Main + interaction effects are estimated as random effects to solve the identifiability problem.
  * Model: Xijk|mu,aj,bk,gjk,Ijk tau~N(mu+aj+bk+Ijk*gjk,tau^-1)
  */

object FPVariableSelection {

  def variableSelection(noOfIter: Int, thin: Int, N: Int, SumObs: Double, structure: DVStructure, alphaPriorMean: Double,  betaPriorMean: Double, mu0: Double, tau0: Double, a: Double, b: Double, thetaPriorMean: Double, aPrior: Double, bPrior: Double, p: Double) = {

    val sampleNo = noOfIter / thin + 1 // Number of samples created from the MCMC
    val alphaLevels= structure.nj
    val betaLevels = structure.nk
    val njk = alphaLevels * betaLevels // Number of levels of interactions

    val curCount = Array(0.0)

    //Define case classes

    case class alphaCoefsCC(acoefs: DenseVector[Double])
    case class betaCoefsCC(bcoefs: DenseVector[Double])
    case class thetaCoefsCC(thcoefs: DenseMatrix[Double])
    case class indicatorsCC(indics: DenseMatrix[Double])
    case class finalCoefsCC(finalCoefs: DenseMatrix[Double])
    case class muTauCC(mt: DenseVector[Double])
    case class tausCC(tauabth: DenseVector[Double])
    case class fullStateCC(acoefs: alphaCoefsCC, bcoefs: betaCoefsCC, thcoefs: thetaCoefsCC, indics: indicatorsCC,finalCoefs: finalCoefsCC, mt: muTauCC, tauabth: tausCC)
    case class fullStateListCC(skta: List[fullStateCC])

    //helper function for mu tau
    def nextmutau(curAlpha: DenseVector[Double], curBeta: DenseVector[Double], curTheta: DenseMatrix[Double], curIndics: DenseMatrix[Double], curmt: muTauCC, curtaus: DenseVector[Double]): muTauCC= {
      val prevtau=curmt.mt(1)
      val prevmu = curmt.mt(0)
      val varMu = 1.0 / (tau0 + N * prevtau) //the variance for mu
      val meanMu = (mu0 * tau0 + prevtau * (SumObs - sumAllMainInterEff(structure, curAlpha, curBeta, alphaLevels, betaLevels, curTheta, curIndics))) * varMu
      val newmu = breeze.stats.distributions.Gaussian(meanMu, sqrt(varMu)).draw()
      val newtau = breeze.stats.distributions.Gamma(a + N / 2.0, 1.0 / (b + 0.5 * YminusMuAndEffects(structure, prevmu, curAlpha, curBeta, curTheta, curIndics))).draw() //  !!!!TO SAMPLE FROM THE GAMMA DISTRIBUTION IN BREEZE THE β IS 1/β
      muTauCC(DenseVector(newmu,newtau))
    }


    //helper function for taus
    def nexttaus(curAlpha: DenseVector[Double] , curBeta: DenseVector[Double], curTheta: DenseMatrix[Double], curIndics: DenseMatrix[Double], curmt: DenseVector[Double], curtaus: tausCC): tausCC= {
      var sumaj = 0.0
      var sumbk = 0.0
      var sumThetajk = 0.0

      for (j <- 0 until alphaLevels) {
        sumaj = sumaj + pow((curAlpha(j) - alphaPriorMean), 2) // Sum used in sampling from Gamma distribution for the precision of alpha
      }

      for (k <- 0 until betaLevels) {
        sumbk = sumbk + pow((curBeta(k) - betaPriorMean), 2) // Sum used in sampling from Gamma distribution for the precision of beta
      }

      for (j <- 0 until alphaLevels) {
        for (k <- 0 until betaLevels) {
          sumThetajk = sumThetajk + pow((curTheta(j, k) - thetaPriorMean), 2) // Sum used in sampling from Gamma distribution for the precision of theta/interacions
        }
      }

      val newtauAlpha = breeze.stats.distributions.Gamma(aPrior + alphaLevels / 2.0, 1.0 / (bPrior + 0.5 * sumaj)).draw() //sample the precision of alpha from gamma
      val newtauBeta = breeze.stats.distributions.Gamma(aPrior + betaLevels / 2.0, 1.0 / (bPrior + 0.5 * sumbk)).draw() // sample the precision of beta from gamma
      val newtauTheta = breeze.stats.distributions.Gamma(aPrior + njk / 2.0, 1.0 / (bPrior + 0.5 * sumThetajk)).draw() // sample the precision of the interactions gamma from gamma Distribition

      tausCC(DenseVector(newtauAlpha, newtauBeta, newtauTheta))
    }

    //helper function for alpha coeffs
    def nextAlphaCoefs(alphaCoefsState: alphaCoefsCC, curBeta: DenseVector[Double], curTheta: DenseMatrix[Double], curIndics: DenseMatrix[Double], curmt: DenseVector[Double], curtaus: DenseVector[Double]): alphaCoefsCC={
      val curAlphaEstim = (DenseVector.zeros[Double](alphaLevels))
      for (j <- 0 until alphaLevels) {
        val SXalphaj = structure.calcAlphaSum(j) // the sum of the observations that have alpha==j
        val Nj = structure.calcAlphaLength(j) // the number of the observations that have alpha==j
        val SumBeta = sumBetaEffGivenAlpha(structure, j, curBeta) //the sum of the beta effects given alpha
        val SinterAlpha = sumInterEffGivenAlpha(structure, j, curTheta, curIndics) //the sum of the gamma/interaction effects given alpha
        val varPalpha = 1.0 / (curtaus(0) + curmt(1) * Nj) //the variance for alphaj
        val meanPalpha = (alphaPriorMean * curtaus(0) + curmt(1) * (SXalphaj - Nj * curmt(0) - SumBeta - SinterAlpha)) * varPalpha //the mean for alphaj
        curAlphaEstim(j) = breeze.stats.distributions.Gaussian(meanPalpha, sqrt(varPalpha)).draw()
      }
      alphaCoefsCC(curAlphaEstim)
    }

    //helper function for beta coeffs
    def nextBetaCoefs(betaCoefsState: betaCoefsCC, curAlpha: DenseVector[Double], curTheta: DenseMatrix[Double], curIndics: DenseMatrix[Double], curmt: DenseVector[Double], curtaus: DenseVector[Double]): betaCoefsCC={
      val curBetaEstim = (DenseVector.zeros[Double](betaLevels))
      for (k <- 0 until betaLevels) {
        val SXbetak = structure.calcBetaSum(k) // the sum of the observations that have beta==k
        val Nk = structure.calcBetaLength(k) // the number of the observations that have beta==k
        val SumAlpha = sumAlphaGivenBeta(structure, k, curAlpha)//the sum of the alpha effects given beta
        val SinterBeta = sumInterEffGivenBeta(structure, k, curTheta, curIndics) //the sum of the gamma/interaction effects given beta
        val varPbeta = 1.0 / (curtaus(1) + curmt(1) * Nk) //the variance for betak
        val meanPbeta = (betaPriorMean * curtaus(1) + curmt(1) * (SXbetak - Nk * curmt(0) - SumAlpha - SinterBeta)) * varPbeta //the mean for betak
        curBetaEstim(k) = breeze.stats.distributions.Gaussian(meanPbeta, sqrt(varPbeta)).draw()
      }
      betaCoefsCC(curBetaEstim)
    }


    //Helper function for indicators and interactions
    def nextIndicsInters(curAlpha: DenseVector[Double], curBeta: DenseVector[Double], curTheta: thetaCoefsCC, curIndics: indicatorsCC,  curmt: DenseVector[Double], curtaus: DenseVector[Double]): (thetaCoefsCC, indicatorsCC,finalCoefsCC )= {
      // Update Indicators and Interaction terms
      val curIndicsEstim = (DenseMatrix.zeros[Double](alphaLevels,betaLevels))
      val curThetaEstim = (DenseMatrix.zeros[Double](alphaLevels,betaLevels))
      var count = 0.0
      for (j <- 0 until alphaLevels) {
        for (k <- 0 until betaLevels) {
          val Njk = structure.getDVList(j, k).length // the number of the observations that have alpha==j and beta==k
          val SXjk = structure.getDVList(j, k).sum // the sum of the observations that have alpha==j and beta==k

          var u = breeze.stats.distributions.Uniform(0, 1).draw()

          //log-sum-exp trick
          val logInitExp = curmt(1) * curTheta.thcoefs(j, k) * (SXjk - Njk * (curmt(0) + curAlpha(j) + curBeta(k) + 0.5 * curTheta.thcoefs(j, k)))
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
            val varPInter = 1.0 / (curtaus(2) + curmt(1) * Njk) //the variance for gammajk
            val meanPInter = (thetaPriorMean * curtaus(2) + curmt(1) * (SXjk - Njk * (curmt(0) + curAlpha(j) + curBeta(k)))) * varPInter
            curThetaEstim(j, k) = breeze.stats.distributions.Gaussian(meanPInter, sqrt(varPInter)).draw()
          }
          else {
            //Update indicator and current interactions if indicator = 0.0
            curIndicsEstim(j, k) = 0.0
            curThetaEstim(j, k) = breeze.stats.distributions.Gaussian(thetaPriorMean, sqrt(1 / curtaus(2))).draw() // sample from the prior of interactions
          }
        }
      }
      curCount(0)= count

      (thetaCoefsCC(curThetaEstim), indicatorsCC(curIndicsEstim), finalCoefsCC(curThetaEstim*:*curIndicsEstim))
    }
    val initmt = muTauCC(DenseVector[Double](0.0,1.0))
    val inittaus = tausCC(DenseVector[Double](1.0,1.0, 1.0))
    val initAlphaCoefs = alphaCoefsCC(DenseVector.zeros[Double](alphaLevels))
    val initBetaCoefs = betaCoefsCC(DenseVector.zeros[Double](betaLevels))



    def calculateNewState( n:Int, fstate:fullStateCC): fullStateCC={

      if (n==0) fstate
      else{
        println(fstate.acoefs.acoefs)
        val latestmt = nextmutau(fstate.acoefs.acoefs, fstate.bcoefs.bcoefs, fstate.thcoefs.thcoefs, fstate.indics.indics, fstate.mt, fstate.tauabth.tauabth)
        val latesttaus = nexttaus(fstate.acoefs.acoefs, fstate.bcoefs.bcoefs, fstate.thcoefs.thcoefs, fstate.indics.indics, latestmt.mt, fstate.tauabth)
        val latestalphas = nextAlphaCoefs(fstate.acoefs, fstate.bcoefs.bcoefs, fstate.thcoefs.thcoefs, fstate.indics.indics, latestmt.mt, latesttaus.tauabth)
        val latestbetas = nextBetaCoefs(fstate.bcoefs, latestalphas.acoefs, fstate.thcoefs.thcoefs, fstate.indics.indics, latestmt.mt, latesttaus.tauabth)
        val (latestthetas, latestIndics, latestfinals) = nextIndicsInters(latestalphas.acoefs, latestbetas.bcoefs, fstate.thcoefs, fstate.indics, latestmt.mt, latesttaus.tauabth )
        val latestFState = fullStateCC(latestalphas, latestbetas, latestthetas,latestIndics, latestfinals, latestmt, latesttaus)
        println(latestFState)
        calculateNewState(noOfIter-1, latestFState)
      }

    }
    calculateNewState(noOfIter, fullStateCC(initAlphaCoefs, initBetaCoefs, thetaCoefsCC(DenseMatrix.zeros[Double](alphaLevels,betaLevels)),indicatorsCC(DenseMatrix.zeros[Double](alphaLevels,betaLevels)), finalCoefsCC(DenseMatrix.zeros[Double](alphaLevels,betaLevels)), initmt, inittaus) )

  }


/**
  * Add all the beta effects for a given alpha.
  */
def sumBetaEffGivenAlpha(structure: DVStructure, alphaIndex: Int, betaEff: DenseVector[Double]): Double = {
  val betaLevels = structure.nk
  var sum = 0.0
  for (ibeta <- 0 until betaLevels) {
    sum += (structure.getDVList(alphaIndex,ibeta).length)*betaEff(ibeta)
  }
  sum
}

/**
  * Add all the alpha effects for a given beta.
  */
def sumAlphaGivenBeta(structure: DVStructure, betaIndex: Int, alphaEff: DenseVector[Double]): Double = {
  val alphaLevels = structure.nj
  var sum = 0.0
  for (ialpha <- 0 until alphaLevels) {
    sum += (structure.getDVList(ialpha,betaIndex).length)*alphaEff(ialpha)
  }
  sum
}
/**
  * Calculate the sum of all the alpha and all the beta effects for all the observations.
  */
def sumAllMainInterEff(structure: DVStructure, alphaEff: DenseVector[Double], betaEff: DenseVector[Double], nj: Int, nk: Int, interEff: DenseMatrix[Double], indics: DenseMatrix[Double]): Double = {
  var sumBeta = 0.0
  var sumAlpha = 0.0
  var sumInter = 0.0
  // For alpha effects
  for (i <- 0 until nk) {
    //through all beta effects
    sumAlpha += sumAlphaGivenBeta(structure, i, alphaEff)
  }

  // For beta effects
  for (i <- 0 until nj) {
    //through all alpha effects
    sumBeta += sumBetaEffGivenAlpha(structure, i, betaEff)
  }

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
  var sum = 0.0
  for (ibeta <- 0 until betaLevels) {
    sum += structure.getDVList(alphaIndex,ibeta).length * indics(alphaIndex, ibeta) * interEff(alphaIndex, ibeta)
  }
  sum
}

/**
  * Add all the interaction effects for a given beta.
  */
def sumInterEffGivenBeta(structure: DVStructure, betaIndex: Int, interEff: DenseMatrix[Double], indics: DenseMatrix[Double]): Double = {
  val alphaLevels = structure.nj
  var sum = 0.0
  for (ialpha <- 0 until alphaLevels) {
    sum += structure.getDVList(ialpha,betaIndex).length * indics(ialpha, betaIndex) * interEff(ialpha, betaIndex)
  }
  sum
}

def sqr(x:DenseVector[Double]) = x*x

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
  val noOfIters = 100000
  val thin = 100
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

  val try1 =
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
}
//  val mtEstim = mean(muTau_est(*, ::)).t
//  val tausEstim = mean(taus_est(*, ::)).t
//  val alphaEstim = mean(alpha_estInter(*, ::)).t
//  val betaEstim = mean(beta_estInter(*, ::)).t
//  val thetaEstim = mean(theta_est(*, ::)).t
//  val indicsEstim = mean(indics_est(*, ::)).t
//  val interactionCoefs = mean(interacs_est(*, ::)).t
//
//  println("Results: ")
//  println("iterations: " + noOfIters)
//  println("thin: " + thin)
//  println("mu and tau:" + mtEstim)
//  println("tauAlpha, tauBeta and tauTheta: " + tausEstim)
//  println("Estimates for alpha: " + alphaEstim)
//  println("Estimates for beta: " + betaEstim)
//  println("Estimates for thetas: " + thetaEstim)
//  println("Estimates for indicators: " + indicsEstim)
//  println("Estimates for interaction coefficients: " + interactionCoefs)
//
//  // Save the results to a csv file
//  val includedM= included_est.toDenseMatrix
//  val mergedMatrix = DenseMatrix.vertcat(muTau_est, taus_est, alpha_estInter, beta_estInter, interacs_est, indics_est,includedM)
//  val outputFIle = new File("/home/antonia/ResultsFromCloud/Report/100919_15x20/Scala/1M/scalaRes1MAfterStr.csv")
//  breeze.linalg.csvwrite(outputFIle, mergedMatrix, separator = ',')


  //    //Plot Results
  //    val resMatR = DenseMatrix.horzcat(muTau_est, alpha_estInter, beta_estInter, theta_est, indics_est)
  //    val plotR = new PlotResults
  //    plotR.plotResults(resMatR, List("Mean", "tau", "alpha 1", "alpha 2", "alpha 3", "beta 1", "beta 2", "beta 3", "beta 4", "a1-b1", "a1-b2", "a1-b3", "a4-b4", "a2-b1", "a2-b2", "a2-b3", "a2-b4", "a3-b1", "a3-b2", "a3-b3", "a3-b4", "Indic a1-b1", "Indic a1-b2", "Indic a1-b3", "Indic a4-b4", "Indic a2-b1", "Indic a2-b2", "Indic a2-b3", "Indic a2-b4", "Indic a3-b1", "Indic a3-b2", "Indic a3-b3", "Indic a3-b4"))


}