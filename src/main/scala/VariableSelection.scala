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

object VariableSelection {

  def variableSelection(noOfIter: Int, thin: Int, N: Int, SumObs: Double, structure: DVStructure, alphaPriorMean: Double,  betaPriorMean: Double, mu0: Double, tau0: Double, a: Double, b: Double, thetaPriorMean: Double, aPrior: Double, bPrior: Double, p: Double) = {

    val sampleNo = noOfIter / thin + 1 // Number of samples created from the MCMC
    val alphaLevels= structure.nj
    val betaLevels = structure.nk
    val njk = alphaLevels * betaLevels // Number of levels of interactions

    val mat_mt = new DenseMatrix[Double](2, sampleNo) // To store mu and tau
    val mat_taus = new DenseMatrix[Double](3, sampleNo) // To store tau_alpha, tau_beta, tau_theta
    val alphaCoefs = new DenseMatrix[Double](alphaLevels, sampleNo) // To store the coefficients for variable alpha
    val betaCoefs = new DenseMatrix[Double](betaLevels, sampleNo) // To store the coefficients for variable beta
    val thetaCoefs = new DenseMatrix[Double](njk, sampleNo) // To store the coefficients for the interactions, variable theta
    val indicators = new DenseMatrix[Double](njk, sampleNo) // To store the indicator variables I
    val finalCoefs= new DenseMatrix[Double](njk, sampleNo) // To store the product of the indicator variables and the estimated coefficient
    val includedInters = new DenseVector[Double](sampleNo) // To store the number of the included interactions at each iteration

    val curAlpha = DenseVector.zeros[Double](alphaLevels) // Current values of the coefficients for every iteration. Initialised with 0.
    val curBeta = DenseVector.zeros[Double](betaLevels)
    val curTheta = DenseMatrix.zeros[Double](alphaLevels, betaLevels)
    val curIndics = DenseMatrix.zeros[Double](alphaLevels, betaLevels)
    val curCount = Array(0.0)

    var mu = 0.0 // Initialise the sampler at mu=0
    var tau = 1.0 // Initialise the sampler at tau=1.0

    var tauAlpha = 1.0 // Initialise the precision for alpha (this is actually tau_alpha)
    var tauBeta = 1.0 // Initialise the precision for beta (this is actually tau_beta)
    var tauTheta = 1.0 // Initialise the precision for theta (this is actually tau_theta)

    var sumaj = 0.0 //Used for sampling the precision of alpha from the FCD (initialised at 0)
    var sumbk = 0.0 //Used for sampling the precision of alpha from the FCD (initialised at 0)
    var sumThetajk = 0.0 //Used for sampling the precision of alpha from the FCD (initialised at 0)

    var ind = 0 //index used for thinning

    //Start the Gibbs sampler from 0s. Initialise with 0 in the first row the matrices where we will store the sampled values.
    mat_mt(::,0) := DenseVector(mu, tau)
    mat_taus(::,0) := DenseVector(tauAlpha, tauBeta, tauTheta)
    alphaCoefs(::,0) := curAlpha
    betaCoefs(::,0) := curBeta
    thetaCoefs(::,0) := curTheta.toDenseVector
    indicators(::,0) := curIndics.toDenseVector
    finalCoefs(::,0) := thetaCoefs(::,0) *:* indicators(::,0)
    includedInters(0) = curCount(0)

    for (i <- 1 until noOfIter) {
      println(i)
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

      tauAlpha = breeze.stats.distributions.Gamma(aPrior + alphaLevels / 2.0, 1.0 / (bPrior + 0.5 * sumaj)).draw() //sample the precision of alpha from gamma
      tauBeta = breeze.stats.distributions.Gamma(aPrior + betaLevels / 2.0, 1.0 / (bPrior + 0.5 * sumbk)).draw() // sample the precision of beta from gamma
      tauTheta = breeze.stats.distributions.Gamma(aPrior + njk / 2.0, 1.0 / (bPrior + 0.5 * sumThetajk)).draw() // sample the precision of the interactions gamma from gamma Distribition

      //Update mu and tau
      val varMu = 1.0 / (tau0 + N * tau) //the variance for mu
      val meanMu = (mu0 * tau0 + tau * (SumObs - sumAllMainInterEff(structure, curAlpha, curBeta, alphaLevels, betaLevels, curTheta, curIndics))) * varMu
      mu = breeze.stats.distributions.Gaussian(meanMu, sqrt(varMu)).draw()
      tau = breeze.stats.distributions.Gamma(a + N / 2.0, 1.0 / (b + 0.5 * YminusMuAndEffects(structure, mu, curAlpha, curBeta, curTheta, curIndics))).draw() //  !!!!TO SAMPLE FROM THE GAMMA DISTRIBUTION IN BREEZE THE β IS 1/β

      // Update alphaj
      for (j <- 0 until alphaLevels) {
        val SXalphaj = structure.calcAlphaSum(j) // the sum of the observations that have alpha==j
        val Nj = structure.calcAlphaLength(j) // the number of the observations that have alpha==j
        val SumBeta = sumBetaEffGivenAlpha(structure, j, curBeta) //the sum of the beta effects given alpha
        val SinterAlpha = sumInterEffGivenAlpha(structure, j, curTheta, curIndics) //the sum of the gamma/interaction effects given alpha
        val varPalpha = 1.0 / (tauAlpha + tau * Nj) //the variance for alphaj
        val meanPalpha = (alphaPriorMean * tauAlpha + tau * (SXalphaj - Nj * mu - SumBeta - SinterAlpha)) * varPalpha //the mean for alphaj
        curAlpha(j) = breeze.stats.distributions.Gaussian(meanPalpha, sqrt(varPalpha)).draw()
      }

      // Update betak
      for (k <- 0 until betaLevels) {
        val SXbetak = structure.calcBetaSum(k) // the sum of the observations that have beta==k
        val Nk = structure.calcBetaLength(k) // the number of the observations that have beta==k
        val SumAlpha = sumAlphaGivenBeta(structure, k, curAlpha)//the sum of the alpha effects given beta
        val SinterBeta = sumInterEffGivenBeta(structure, k, curTheta, curIndics) //the sum of the gamma/interaction effects given beta
        val varPbeta = 1.0 / (tauBeta + tau * Nk) //the variance for betak
        val meanPbeta = (betaPriorMean * tauBeta + tau * (SXbetak - Nk * mu - SumAlpha - SinterBeta)) * varPbeta //the mean for betak
        curBeta(k) = breeze.stats.distributions.Gaussian(meanPbeta, sqrt(varPbeta)).draw()
      }

      // Update Indicators and Interaction terms
      var count = 0.0
      for (j <- 0 until alphaLevels) {
        for (k <- 0 until betaLevels) {
          val Njk = structure.getDVList(j, k).length // the number of the observations that have alpha==j and beta==k
          val SXjk = structure.getDVList(j, k).sum // the sum of the observations that have alpha==j and beta==k

          var u = breeze.stats.distributions.Uniform(0, 1).draw()

          //log-sum-exp trick
          val logInitExp = tau * curTheta(j, k) * (SXjk - Njk * (mu + curAlpha(j) + curBeta(k) + 0.5 * curTheta(j, k)))
          val logProb0 = log(1.0 - p) //The log of the probability I=0
          val logProb1 = log(p) + logInitExp //The log of the probability I=1
          val maxProb = max(logProb0, logProb1) //Find the max of the two probabilities
          val scaledProb0 = exp(logProb0 - maxProb) //Scaled by subtracting the max value and exponentiating
          val scaledProb1 = exp(logProb1 - maxProb) //Scaled by subtracting the max value and exponentiating
          var newProb0 = scaledProb0 / (scaledProb0 + scaledProb1) //Normalised
          val newProb1 = scaledProb1 / (scaledProb0 + scaledProb1) //Normalised

          if (newProb0 < u) {
            //prob0: Probability for when the indicator = 0, so if prob0 < u => indicator = 1
            curIndics(j, k) = 1.0
            count += 1.0
            val varPInter = 1.0 / (tauTheta + tau * Njk) //the variance for gammajk
            val meanPInter = (thetaPriorMean * tauTheta + tau * (SXjk - Njk * (mu + curAlpha(j) + curBeta(k)))) * varPInter
            curTheta(j, k) = breeze.stats.distributions.Gaussian(meanPInter, sqrt(varPInter)).draw()
          }
          else {
            //Update indicator and current interactions if indicator = 0.0
            curIndics(j, k) = 0.0
            curTheta(j, k) = breeze.stats.distributions.Gaussian(thetaPriorMean, sqrt(1 / tauTheta)).draw() // sample from the prior of interactions
          }
        }
      }
      curCount(0)= count

      // Thinning
      if ((i % thin).equals(0)) {
        ind = ind + 1
        mat_mt(::, ind) := DenseVector(mu, tau)
        mat_taus(::, ind) := DenseVector(tauAlpha, tauBeta, tauTheta)
        alphaCoefs(::, ind) := curAlpha
        betaCoefs(::, ind) := curBeta
        thetaCoefs(::, ind) := curTheta.toDenseVector
        indicators(::, ind) := curIndics.toDenseVector
        finalCoefs(::, ind):= thetaCoefs(::, ind) *:* indicators(::, ind)
        includedInters(ind) = curCount(0)
      }

      sumaj = 0.0
      sumbk = 0.0
      sumThetajk = 0.0
    }
    finalCoefs:= thetaCoefs *:* indicators
    (mat_mt, mat_taus, alphaCoefs, betaCoefs, thetaCoefs, indicators, finalCoefs, includedInters)
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
    val alpha = data(::, 1).map(_.toInt).map(x=>x-1)
    val beta = data(::, 2).map(_.toInt).map(x=>x-1)
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

    val (muTau_est, taus_est, alpha_estInter, beta_estInter, theta_est, indics_est, interacs_est, included_est) =
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

    val mtEstim = mean(muTau_est(*, ::)).t
    val tausEstim = mean(taus_est(*, ::)).t
    val alphaEstim = mean(alpha_estInter(*, ::)).t
    val betaEstim = mean(beta_estInter(*, ::)).t
    val thetaEstim = mean(theta_est(*, ::)).t
    val indicsEstim = mean(indics_est(*, ::)).t
    val interactionCoefs = mean(interacs_est(*, ::)).t

    println("Results: ")
    println("iterations: " + noOfIters)
    println("thin: " + thin)
    println("mu and tau:" + mtEstim)
    println("tauAlpha, tauBeta and tauTheta: " + tausEstim)
    println("Estimates for alpha: " + alphaEstim)
    println("Estimates for beta: " + betaEstim)
    println("Estimates for thetas: " + thetaEstim)
    println("Estimates for indicators: " + indicsEstim)
    println("Estimates for interaction coefficients: " + interactionCoefs)

    // Save the results to a csv file
    val includedM= included_est.toDenseMatrix
    val mergedMatrix = DenseMatrix.vertcat(muTau_est, taus_est, alpha_estInter, beta_estInter, interacs_est, indics_est,includedM)
    val outputFIle = new File("/home/antonia/ResultsFromCloud/Report/100919_15x20/Scala/1M/scalaRes1MAfterStr.csv")
    breeze.linalg.csvwrite(outputFIle, mergedMatrix, separator = ',')


    //    //Plot Results
    //    val resMatR = DenseMatrix.horzcat(muTau_est, alpha_estInter, beta_estInter, theta_est, indics_est)
    //    val plotR = new PlotResults
    //    plotR.plotResults(resMatR, List("Mean", "tau", "alpha 1", "alpha 2", "alpha 3", "beta 1", "beta 2", "beta 3", "beta 4", "a1-b1", "a1-b2", "a1-b3", "a4-b4", "a2-b1", "a2-b2", "a2-b3", "a2-b4", "a3-b1", "a3-b2", "a3-b3", "a3-b4", "Indic a1-b1", "Indic a1-b2", "Indic a1-b3", "Indic a4-b4", "Indic a2-b1", "Indic a2-b2", "Indic a2-b3", "Indic a2-b4", "Indic a3-b1", "Indic a3-b2", "Indic a3-b3", "Indic a3-b4"))
  }

}