package mcmc.gibbs

import java.io.{File, FileWriter, PrintWriter}
import breeze.linalg.{DenseMatrix, DenseVector, max}
import breeze.numerics.{exp, log, pow, sqrt}
import structure.DVStructure
import scala.math.Pi

/**
  * Variable selection with Gibbs sampler. Implementation for symmetric main effects and asymmetric interactions.
  * Model: X_ijk | mu,a_j,b_k ,I_jk,theta_jk,tau  ~ N(mu + z_j + z_k + I_jk * theta_jk , τ^−1 )
  * Using gamma priors for taua, taub, tauTheta, Bermouli for the variable selection indicator I_jk and Normal for the main effects z and the effect size theta_jk
  * Symmetric main effects: zs come from the same distribution
  * Asymmetric Interactions: I_jk * theta_jk != I_kj * theta_kj
  **/
class SymmetricMain extends VariableSelection {

  override def variableSelection(info: InitialInfo) = {
    // Initialise case class objects
    val initmt = DenseVector[Double](0.0,1.0)
    val inittaus = DenseVector[Double](1.0,1.0)
    val initAlphaCoefs = DenseVector.zeros[Double](info.alphaLevels) //Not used in SymmetricMain implementation
    val initBetaCoefs = DenseVector.zeros[Double](info.betaLevels) //Not used in SymmetricMain implementation
    val initZetaCoefs = DenseVector.zeros[Double](info.zetaLevels)
    val initThetas = DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels) //zetaLevels in SymmetricMain implementation
    val initIndics = DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels) //zetaLevels in SymmetricMain implementation
    val initFinals = DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels) //zetaLevels in SymmetricMain implementation
    val loglik = 0.0
    val incProb = 0.0

    val fullStateInit = FullState(initAlphaCoefs, initBetaCoefs, initZetaCoefs, initThetas, initIndics, initFinals, initmt, inittaus, loglik, incProb)
    calculateAllStates(info.noOfIter, info, fullStateInit)
  }

  /**
    * Function for updating the inclusion probability p
    */
  override def nextInclusionProb(oldfullState: FullState, info: InitialInfo): FullState = {
    val njk = info.structure.sizeOfStructure()
    val SumIjk = breeze.linalg.sum(oldfullState.indics)
    val alphaP = info.pap + SumIjk
    val betaP = info.pbp + njk - SumIjk
    val newp = breeze.stats.distributions.Beta.distribution(alphaP, betaP).draw()
    oldfullState.copy(inclProb = newp)
  }

  /**
    * Function for updating mu and tau
    */
  override def nextmutau(oldfullState: FullState, info: InitialInfo): FullState= {
    val varMu = 1.0 / (info.tau0 + info.N * oldfullState.mt(1)) //the variance for mu
    val meanMu = (info.mu0 * info.tau0 + oldfullState.mt(1) * (info.SumObs - sumAllMainInterEff(info.structure, oldfullState.zcoefs, info.zetaLevels, oldfullState.thcoefs, oldfullState.indics))) * varMu
    val newmu = breeze.stats.distributions.Gaussian(meanMu, sqrt(varMu)).draw()
    val obsMinusFitted = YminusMuAndEffects(info.structure, newmu, oldfullState.zcoefs, oldfullState.thcoefs, oldfullState.indics)
    //Use the just updated mu to estimate tau
    val newtau = breeze.stats.distributions.Gamma(info.a + info.N / 2.0, 1.0 / (info.b + 0.5 * obsMinusFitted)).draw() //  !!!!TO SAMPLE FROM THE GAMMA DISTRIBUTION IN BREEZE THE β IS 1/β
    val newLogLik = if(info.logLikFlag) {
      0.5 * info.N * (log(newtau) - log(2 * Pi)) - 0.5 * newtau * obsMinusFitted // Caclulate the logLikelihood
    }else{
      0.0
    }
    oldfullState.copy(mt=DenseVector(newmu,newtau), logLik = newLogLik)
  }

  /**
    * Function for updating taus (tauz, tauInt)
    */
  override def nexttaus(oldfullState: FullState, info: InitialInfo):FullState= {
    val njk = info.structure.sizeOfStructure() // Number of levels of interactions

    var sumzj = 0.0

    oldfullState.zcoefs.foreachValue( zcoef => {
      sumzj += pow(zcoef - info.alphaPriorMean, 2)
    })
    sumzj -= (info.zetaLevels - info.zetaLevelsDist) * pow(0 - info.alphaPriorMean, 2) //For the missing effects (if any) added extra in the sum above

    var sumThetajk = 0.0
    oldfullState.thcoefs.foreachValue(thcoef => {
      sumThetajk += pow(thcoef -info.thetaPriorMean, 2) // Sum used in sampling from Gamma distribution for the precision of theta/interacions
    })

    sumThetajk -= (info.zetaLevels * info.zetaLevels - njk) * pow(0 - info.betaPriorMean, 2) //For the missing effects (if any) added extra in the sum above

    val newtauZeta = breeze.stats.distributions.Gamma(info.aPrior + info.zetaLevels / 2.0, 1.0 / (info.bPrior + 0.5 * sumzj)).draw() //sample the precision of alpha from gamma
    val newtauTheta = breeze.stats.distributions.Gamma(info.aPrior + njk / 2.0, 1.0 / (info.bPrior + 0.5 * sumThetajk)).draw() // sample the precision of the interactions gamma from gamma Distribition

    oldfullState.copy(tauabth = DenseVector(newtauZeta, newtauTheta))
  }

  override def nextCoefs(oldfullState: FullState, info: InitialInfo): FullState = {
    nextZetaCoefs(oldfullState, info)
  }

  /**
    * Function for updating zeta coefficients.
    * Each zeta depends on the other zs, for which the latest update needs to be used.
    */
  def nextZetaCoefs(oldfullState: FullState, info: InitialInfo):FullState={

    val curZetaEstim = DenseVector.zeros[Double](info.zetaLevels)
    curZetaEstim:= oldfullState.zcoefs

    info.structure.getAllZetas().foreach( item => { //For each existing zeta
      val j = item
      val SXZetaj = info.structure.calcZetaSum(j) // the sum of the observations that have zeta == j on either side, not both
      val Nj = info.structure.calcZetaLength(j) // the number of the observations that have zeta == j on either side, not both
      val Njj = info.structure.calcDoubleZetaLength(j) // the number of the observations that have zeta == j on both sides
      val SXZetajDouble = info.structure.calcDoubleZetaSum(j) // the sum of the observations that have zeta == j on both sides
      val SumZeta = sumEffectsOfOtherZetas(info.structure, j, curZetaEstim) //the sum of the other zeta effects given zeta, for which the given z is on either side (but not on both sides)
      val SinterZeta = sumInterEffGivenZeta(info.structure, j, oldfullState.thcoefs, oldfullState.indics) //the sum of the gamma/interaction effects given zeta, for which the given z is on either side (but not on both sides)
      val SinterZetaDoubles = sumInterEffDoublesGivenZeta(info.structure, j, oldfullState.thcoefs, oldfullState.indics) //the sum of the gamma/interaction effects given zeta, for which the given z is on both sides
      val varPzeta = 1.0 / (oldfullState.tauabth(0) + oldfullState.mt(1) * Nj + 4 * oldfullState.mt(1) * Njj) //the variance for zetaj
      val meanPzeta = (info.alphaPriorMean * oldfullState.tauabth(0) + oldfullState.mt(1) * (SXZetaj - Nj * oldfullState.mt(0) - SumZeta - SinterZeta + 2 * SXZetajDouble - 2 * Njj * oldfullState.mt(0) - 2 * SinterZetaDoubles )) * varPzeta //the mean for alphaj
      curZetaEstim.update(j, breeze.stats.distributions.Gaussian(meanPzeta, sqrt(varPzeta)).draw())
    })

    oldfullState.copy(zcoefs = curZetaEstim)
  }

  /**
    * Function for updating indicators, interactions and final interaction coefficients
    */
  override def nextIndicsInters(oldfullState: FullState, info: InitialInfo):FullState= {

    val estimations = info.structure.map(item => ((item.a, item.b), {
      val Njk = item.list.length // the number of the observations that have alpha==j and beta==k
      val SXjk = item.list.sum // the sum of the observations that have alpha==j and beta==k

      val u = breeze.stats.distributions.Uniform(0, 1).draw()

      //log-sum-exp trick
      val thcoef = oldfullState.thcoefs(item.a, item.b)
      val logInitExp = oldfullState.mt(1) * thcoef * (SXjk - Njk * (oldfullState.mt(0) + oldfullState.zcoefs(item.a) + oldfullState.zcoefs(item.b) + 0.5 * thcoef))
      val logProb0 = log(1.0 - oldfullState.inclProb) //The log of the probability I=0
      val logProb1 = log(oldfullState.inclProb) + logInitExp //The log of the probability I=1
      val maxProb = max(logProb0, logProb1) //Find the max of the two probabilities
      val scaledProb0 = exp(logProb0 - maxProb) //Scaled by subtracting the max value and exponentiating
      val scaledProb1 = exp(logProb1 - maxProb) //Scaled by subtracting the max value and exponentiating
      val newProb0 = scaledProb0 / (scaledProb0 + scaledProb1) //Normalised
      // val newProb1 = scaledProb1 / (scaledProb0 + scaledProb1) //Normalised

      if (newProb0 < u) {
        //prob0: Probability for when the indicator = 0, so if prob0 < u => indicator = 1
        val indicsEstim = 1.0
        val varPInter = 1.0 / (oldfullState.tauabth(1) + oldfullState.mt(1) * Njk) //the variance for gammajk
        val meanPInter = (info.thetaPriorMean * oldfullState.tauabth(1) + oldfullState.mt(1) * (SXjk - Njk * (oldfullState.mt(0) + oldfullState.zcoefs(item.a) + oldfullState.zcoefs(item.b)))) * varPInter
        val thetaEstim = breeze.stats.distributions.Gaussian(meanPInter, sqrt(varPInter)).draw()
        (indicsEstim, thetaEstim)
      }
      else {
        //Update indicator and current interactions if indicator = 0.0
        val indicsEstim = 0.0
        val thetaEstim = breeze.stats.distributions.Gaussian(info.thetaPriorMean, sqrt(1 / oldfullState.tauabth(1))).draw() // sample from the prior of interactions
        (indicsEstim, thetaEstim)
      }
    }))
    val curIndicsEstim = DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels)
    val curThetaEstim = DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels)

    estimations.seq //make sure we are on single thread to access the collections above without concurrency problem
      .foreach { case (key, value) =>
        curIndicsEstim(key._1, key._2) = value._1
        curThetaEstim(key._1, key._2) = value._2
      }

    oldfullState.copy(thcoefs = curThetaEstim, indics = curIndicsEstim, finalCoefs = curThetaEstim*:*curIndicsEstim)
  }

  /**
    * Add all the zeta effects for all the other zetas for that specific zeta.
    * e.g. updating z1: (1,1),(1,2),(2,1),(1,3),(1,4),(4,1) => Sum the effects for: z2*NoOfObs for that category + z2*NoOfObs for that category + z3*NoOfObs for that category + z4*NoOfObs for that category + z4*NoOfObs for that category
    */
  def sumEffectsOfOtherZetas(structure: DVStructure, zetaIndex: Int, zetaEff: DenseVector[Double]): Double = {
    //returns the element which is not zetaIndex. It doesn't take into account the cases where both sides are zetaIndex because getAllOtherZetasItemsForGivenZ works on a structure that does not involve the (j,j) cases

    //Alternative implementations, more functional but less efficient
    //structure.getAllOtherZetasItemsForGivenZ(zetaIndex).map(elem => elem._2.length * zetaEff(notZeta(elem._1._1, elem._1._2))).sum
    //structure.getAllOtherZetasItemsForGivenZ(zetaIndex).foldLeft(0.0)( (sum, elem) => sum + (elem._2.length * zetaEff(notZeta(elem._1._1, elem._1._2))) )

    var totalsum = 0.0
    structure.getAllOtherZetasItemsForGivenZ(zetaIndex).foreach(item => {
      val otherElement = if (item._1._1 == zetaIndex) item._1._2 else item._1._1
      totalsum += item._2.length * zetaEff(otherElement)
    })
    totalsum
  }

  /**
    * Calculate the sum of all the zeta 1 and all the zeta 2 effects for all the observations.
    */
  def sumAllMainInterEff(structure: DVStructure, zetaEff: DenseVector[Double], nz: Int, interEff: DenseMatrix[Double], indics: DenseMatrix[Double]): Double = {
    structure.map(item => {
      item.list.length * (zetaEff(item.a) + zetaEff(item.b) + indics(item.a, item.b) * interEff(item.a, item.b))
    }).sum
  }

  /**
    * Add all the interaction effects for a given zeta. Adds all the interactions for which zeta is on either side. Includes the doubles bcs getZetasItemsForGivenZ uses a structure that includes everything
    */
  def sumInterEffGivenZeta(structure: DVStructure, zetaIndex: Int, interEff: DenseMatrix[Double], indics: DenseMatrix[Double]): Double = {
    //Alternative implementations more functional, but less efficient. foldLeft faster than map (50x60, 100K, 20.3sec vs 46 sec). var total sum the fastest (16.4 sec)
    //structure.getAllOtherZetasItemsForGivenZ(zetaIndex).map(elem => elem._2.length * indics(elem._1._1, elem._1._2) * interEff(elem._1._1, elem._1._2)).sum
    //structure.getAllOtherZetasItemsForGivenZ(zetaIndex).foldLeft(0.0)( (sum, elem) => sum + (elem._2.length * indics(elem._1._1, elem._1._2) * interEff(elem._1._1, elem._1._2)) )

    var totalsum = 0.0
    structure.getAllOtherZetasItemsForGivenZ(zetaIndex).foreach(item => {
      totalsum += item._2.length * indics(item._1._1, item._1._2) * interEff(item._1._1, item._1._2)
    })
    totalsum
  }

  /**
    * Add all the interaction effects for a given zeta which is double (zeta,zeta)
    */
  def sumInterEffDoublesGivenZeta(structure: DVStructure, zetaIndex: Int, interEff: DenseMatrix[Double], indics: DenseMatrix[Double]): Double = {
    structure.getAllDoubleZetasItemsForGivenZ(zetaIndex).foldLeft(0.0)( (sum, elem) => sum + (elem._2.length * indics(elem._1._1, elem._1._2) * interEff(elem._1._1, elem._1._2)) )
  }

  /**
    * Calculate the Yi-mu-u_eff-n_eff- inter_effe. To be used in estimating tau
    */
  def YminusMuAndEffects(structure:DVStructure, mu: Double, zetaEff: DenseVector[Double], interEff: DenseMatrix[Double], indics: DenseMatrix[Double]): Double = {
      structure.map( item => {
        val a = item.a
        val b = item.b
        item.list.foldLeft(0.0)((sum, x) => sum + scala.math.pow(x - mu - zetaEff(a) - zetaEff(b) - interEff(a, b) * indics(a, b), 2))
      }).sum
  }

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

    pw.append("mu ,tau, tauz, tauInt, logLik, p,")
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
        .append(fullstate.logLik.toString)
        .append(",")
        .append(fullstate.inclProb.toString)
        .append(",")
        .append( fullstate.zcoefs.toArray.map { alpha => alpha.toString }.mkString(",") )
        .append(",")
        .append( fullstate.finalCoefs.toArray.map { theta => theta.toString }.mkString(",") )
        .append(",")
        .append( fullstate.indics.toArray.map { ind => ind.toString }.mkString(",") )
        .append("\n")
    }
    pw.close()
  }

}
