import java.io.File
import breeze.linalg.{*, _}
import breeze.numerics._
import cats._
import breeze.stats.mean

object MainSymmetricNew {
  def variableSelection(noOfIter: Int, thin: Int, N: Int, SumObs: Double, structure: DVStructure, alphaLevels: Int, betaLevels: Int, zetaLevels: Int,
                        alphaPriorMean: Double, betaPriorMean: Double, thetaPriorMean: Double, mu0: Double, tau0: Double,
                        a: Double, b: Double, aPrior: Double, bPrior: Double, p: Double) = {

    val njk = alphaLevels * betaLevels // Number of levels of interactions

    val curCount = Array(0.0)

    //Define case classes
    case class FullState(zcoefs: DenseVector[Double], thcoefs: DenseMatrix[Double], indics: DenseMatrix[Double],finalCoefs: DenseMatrix[Double], mt: DenseVector[Double], tauabth: DenseVector[Double])

    case class FullStateList(fstateL: List[FullState])

    // Update mu and tau
    // helper function for mu tau
    def nextmutau(oldfullState: FullState): FullState= {
      val prevtau = oldfullState.mt(1)
      val prevmu = oldfullState.mt(0)
      val varMu = 1.0 / (tau0 + N * prevtau) //the variance for mu
      val meanMu = (mu0 * tau0 + prevtau * (SumObs - sumAllMainInterEff(structure, oldfullState.zcoefs, zetaLevels, oldfullState.thcoefs, oldfullState.indics))) * varMu
      val newmu = breeze.stats.distributions.Gaussian(meanMu, sqrt(varMu)).draw()
      val newtau = breeze.stats.distributions.Gamma(a + N / 2.0, 1.0 / (b + 0.5 * YminusMuAndEffects(structure, prevmu, oldfullState.zcoefs, oldfullState.thcoefs, oldfullState.indics))).draw() //  !!!!TO SAMPLE FROM THE GAMMA DISTRIBUTION IN BREEZE THE β IS 1/β
      oldfullState.copy(mt=DenseVector(newmu,newtau))
    }

    // Update taus (taua, taub, tauInt)
    // helper function for taus
    def nexttaus(oldfullState: FullState):FullState= {

      //todo: check if acoef non set values create an issue
      var sumzj = 0.0
      oldfullState.zcoefs.foreachValue( zcoef => {
        sumzj += pow(zcoef - alphaPriorMean, 2)
      })

      //todo: check if thcoef non set values create an issue
      var sumThetajk = 0.0
      oldfullState.thcoefs.foreachValue(thcoef => {
        sumThetajk += pow(thcoef - thetaPriorMean, 2) // Sum used in sampling from Gamma distribution for the precision of theta/interacions
      })

      val newtauZeta = breeze.stats.distributions.Gamma(aPrior + zetaLevels / 2.0, 1.0 / (bPrior + 0.5 * sumzj)).draw() //sample the precision of alpha from gamma
      val newtauTheta = breeze.stats.distributions.Gamma(aPrior + njk / 2.0, 1.0 / (bPrior + 0.5 * sumThetajk)).draw() // sample the precision of the interactions gamma from gamma Distribition

      oldfullState.copy(tauabth = DenseVector(newtauZeta, newtauTheta))
    }

    // Update Zeta coefficients
    // helper function for Zeta coeffs
    def nextZetaCoefs(oldfullState: FullState):FullState={
      val curZetaEstim = (DenseVector.zeros[Double](zetaLevels))
      structure.getAllItemsMappedByZ().foreach( item => {
        val j = item._1
        val SXzetaj = structure.calcZetaSum(j) // the sum of the observations that have zeta on either side
        val Nj = structure.calcZetaLength(j) // the sum of the observations that have zeta on either side
        val SumZeta = sumRemZetaEffGivenZeta(structure, j, oldfullState.zcoefs) //the sum of the zeta effects of the rest for a given zeta
        val SinterZeta = sumInterEffGivenZeta(structure, j, oldfullState.thcoefs, oldfullState.indics) //the sum of the gamma/interaction effects given zeta
        val varPzeta = 1.0 / (oldfullState.tauabth(0) + oldfullState.mt(1) * Nj) //the variance for zetaj
        val meanPzeta = (alphaPriorMean * oldfullState.tauabth(0) + oldfullState.mt(1) * (SXzetaj - Nj * oldfullState.mt(0) - SumZeta - SinterZeta)) * varPzeta //the mean for zetaj
        curZetaEstim(j) = breeze.stats.distributions.Gaussian(meanPzeta, sqrt(varPzeta)).draw()
      })
      oldfullState.copy(zcoefs = curZetaEstim)
    }

    // Update indicators, interactions and final interaction coefficients
    //Helper function for indicators and interactions
    def nextIndicsInters(oldfullState: FullState):FullState= {

      val curIndicsEstim = (DenseMatrix.zeros[Double](alphaLevels,betaLevels))
      val curThetaEstim = (DenseMatrix.zeros[Double](alphaLevels,betaLevels))
      var count = 0.0

      structure.foreach( item => {
        val Njk = item.list.length // the number of the observations that have alpha==j and beta==k
        val SXjk = item.list.sum // the sum of the observations that have alpha==j and beta==k

        val u = breeze.stats.distributions.Uniform(0, 1).draw()

        //log-sum-exp trick
        val thcoef = oldfullState.thcoefs(item.a, item.b)
        val logInitExp = oldfullState.mt(1) * thcoef * (SXjk - Njk * (oldfullState.mt(0) + oldfullState.zcoefs(item.a) + oldfullState.zcoefs(item.b) + 0.5 * thcoef))
        val logProb0 = log(1.0 - p) //The log of the probability I=0
        val logProb1 = log(p) + logInitExp //The log of the probability I=1
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
          val meanPInter = (thetaPriorMean * oldfullState.tauabth(1) + oldfullState.mt(1) * (SXjk - Njk * (oldfullState.mt(0) + oldfullState.zcoefs(item.a) + oldfullState.zcoefs(item.b)))) * varPInter
          curThetaEstim(item.a, item.b) = breeze.stats.distributions.Gaussian(meanPInter, sqrt(varPInter)).draw()
        }
        else {
          //Update indicator and current interactions if indicator = 0.0
          curIndicsEstim(item.a,item.b) = 0.0
          curThetaEstim(item.a,item.b) = breeze.stats.distributions.Gaussian(thetaPriorMean, sqrt(1 / oldfullState.tauabth(1))).draw() // sample from the prior of interactions
        }
      })

      curCount(0)= count
      oldfullState.copy(thcoefs = curThetaEstim, indics = curIndicsEstim, finalCoefs = curThetaEstim*:*curIndicsEstim)
    }

    // Initialise case class objects
    val initmt = DenseVector[Double](0.0,1.0)
    val inittaus = DenseVector[Double](1.0,1.0)
    val initZetaCoefs = DenseVector.zeros[Double](zetaLevels)
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
        println(n)
        val latestmt = nextmutau(fstate)
        val latesttaus = nexttaus(latestmt)
        val latestzetas = nextZetaCoefs(latesttaus)
        val latestFullyUpdatedState = nextIndicsInters(latestzetas)
        if((n % thin).equals(0)) {
          calculateNewState(n-1, latestFullyUpdatedState, FullStateList(latestFullyUpdatedState::fstateList.fstateL))
        }
        else calculateNewState(n-1, latestFullyUpdatedState, fstateList)
      }
    }
    calculateNewState(noOfIter, FullState(initZetaCoefs, initThetas, initIndics, initFinals, initmt, inittaus), FullStateList(List(FullState(initZetaCoefs, initThetas, initIndics, initFinals, initmt, inittaus))))
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
    val data = csvread(new File("/home/antonia/ResultsFromCloud/Report/Symmetric/symmetricMain/simulInterSymmetricMain.csv"))
    val sampleSize = data.rows
    val y = data(::, 0)
    val sumObs = y.toArray.sum // Sum of the values of all the observations
    val alpha = data(::, 1).map(_.toInt).map(x => x - 1)
    val beta = data(::, 2).map(_.toInt).map(x => x - 1)
    //    val structure : DVStructure = new DVStructureArrays(y, alpha, beta)
    val structure : DVStructure = new DVStructureMap(y, alpha, beta)
    val alphaLevels = alpha.toArray.distinct.length
    val betaLevels = beta.toArray.distinct.length
    val zetaLevels = max(alphaLevels, betaLevels)

    // Parameters
    val noOfIters = 100000
    val thin = 10
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
          noOfIters, thin, sampleSize, sumObs, structure, alphaLevels, betaLevels, zetaLevels,
          alphaPriorMean, betaPriorMean, interPriorMean, mu0, tau0,
          a, b, aPrior, bPrior, p)
      )

    println("zetas")
    val bcoefficients = statesResults.fstateL.map(f=>f.zcoefs)
    val bcoefMat= DenseMatrix(bcoefficients.map(_.toArray):_*)
    val meanValsBcoef = mean(bcoefMat(::, *))
    println(meanValsBcoef)

    println("mu, tau")
    val mtcoefficients = statesResults.fstateL.map(f=>f.mt)
    val mtcoefMat= DenseMatrix(mtcoefficients.map(_.toArray):_*)
    val meanValsmtcoef = mean(mtcoefMat(::, *))
    println(meanValsmtcoef)

    println("taus")
    val tauscoefficients = statesResults.fstateL.map(f=>f.tauabth)
    val tauscoefMat= DenseMatrix(tauscoefficients.map(_.toArray):_*)
    //val outputFile = new File("/home/antonia/Desktop/tausTry.csv")
    //breeze.linalg.csvwrite(outputFile, tauscoefMat, separator = ',')
    val meanValstauscoef = mean(tauscoefMat(::, *))
    println(meanValstauscoef)

    println("thetas")
    val thetascoefficients = statesResults.fstateL.map(f=>f.thcoefs.toDenseVector)
    val thetascoefMat= DenseMatrix(thetascoefficients.map(_.toArray):_*)
    val meanValsthetascoef = mean(thetascoefMat(::, *))
    println(meanValsthetascoef)

    println("indicators")
    val indicscoefficients = statesResults.fstateL.map(f=>f.indics.toDenseVector)
    val indicscoefMat= DenseMatrix(indicscoefficients.map(_.toArray):_*)
    val meanValsindicscoef = mean(indicscoefMat(::, *))
    println(meanValsindicscoef)

    println("finalCoefs")
    val finalcoefficients = statesResults.fstateL.map(f=>f.finalCoefs.toDenseVector)
    val finalcoefMat= DenseMatrix(finalcoefficients.map(_.toArray):_*)
    val meanValsfinalcoef = mean(finalcoefMat(::, *))
    println(meanValsfinalcoef)
  }

}
