package mcmc.gibbs

import breeze.linalg.{DenseMatrix, DenseVector, max, upperTriangular}
import breeze.numerics.{exp, log, pow, sqrt}

class SymmetricBoth extends SymmetricMain3 {

  override def variableSelection(info: InitialInfo): FullStateList = {
    // Initialise case class objects
    val initmt = DenseVector[Double](0.0,1.0)
    val inittaus = DenseVector[Double](1.0,1.0)
    val initAlphaCoefs = DenseVector.zeros[Double](info.alphaLevels)
    val initBetaCoefs = DenseVector.zeros[Double](info.betaLevels)
    val initZetaCoefs = DenseVector.zeros[Double](info.zetaLevels)
    val initThetas = DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels)
    val initIndics = DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels)
    val initFinals = DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels)

    calculateNewState(info.noOfIter, info, FullState(initAlphaCoefs, initBetaCoefs, initZetaCoefs, initThetas, initIndics, initFinals, initmt, inittaus), FullStateList(List(FullState(initAlphaCoefs, initBetaCoefs, initZetaCoefs, initThetas, initIndics, initFinals, initmt, inittaus))))
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
    upperTriangular(oldfullState.thcoefs).foreachValue(thcoef => {
      sumThetajk += pow(thcoef - info.thetaPriorMean, 2) // Sum used in sampling from Gamma distribution for the precision of theta/interacions
    })
    sumThetajk += (info.noOfInters - info.sizeOfDouble)*pow(- info.thetaPriorMean, 2)
    val njk = info.noOfInters // Number of levels of interactions
    val newtauZeta = breeze.stats.distributions.Gamma(info.aPrior + info.zetaLevels / 2.0, 1.0 / (info.bPrior + 0.5 * sumzj)).draw() //sample the precision of alpha from gamma
    val newtauTheta = breeze.stats.distributions.Gamma(info.aPrior + njk / 2.0, 1.0 / (info.bPrior + 0.5 * sumThetajk)).draw() // sample the precision of the interactions gamma from gamma Distribition

    oldfullState.copy(tauabth = DenseVector(newtauZeta, newtauTheta))
    //oldfullState.copy(tauabth = DenseVector(0.385, 0.265))
  }

  // Update indicators, interactions and final interaction coefficients
  //Helper function for indicators and interactions
  override def nextIndicsInters(oldfullState: FullState, info: InitialInfo): FullState= {

    val curIndicsEstim = (DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels))
    val curThetaEstim = (DenseMatrix.zeros[Double](info.zetaLevels, info.zetaLevels))
    var count = 0.0

    info.structureSorted.foreach( item => {
      val j= item.a
      val k = item.b
      // toDo: Njkkj is the same as item.list.length from structureSorted so see maybe it is not necessary to have calcAlphaBetaLength
      // the number of the observations that have alpha==j and beta==k and alpha==k and beta==j
      val Njkkj =
      if(j==k) {info.structure.calcAlphaBetaLength(j,k)
      }else{
        info.structure.calcAlphaBetaLength(j,k) + info.structure.calcAlphaBetaLength(k,j)
      }

      // the number of the observations that have alpha==j and beta==k and alpha==k and beta==j
      // the sum of the observations that have alpha==j and beta==k and alpha==k and beta==j
      val SXjkkj =
      if(j==k) {
        info.structure.calcAlphaBetaSum(j,k)
      } else{
        info.structure.calcAlphaBetaSum(j,k) + info.structure.calcAlphaBetaSum(k,j)
      }

      val u = breeze.stats.distributions.Uniform(0, 1).draw()

      //log-sum-exp trick
      val thcoef = oldfullState.thcoefs(item.a, item.b) // This is the same for (item.a, item.b) and (item.b, item.a) so it does not matter which one we use
      val NoOfajForbk = info.structure.calcAlphaBetaLength(j,k) //No of observations for which a==j and b==k
      val NoOfakForbj = info.structure.calcAlphaBetaLength(k,j) //No of observations for which a==k and b==j

      def returnIfExists(dv: DenseVector[Double], ind: Int)={
        if(ind < dv.length) dv(ind)
        else 0.0
      }

      val SigmaTheta =
        if(j==k){
          SXjkkj - Njkkj * oldfullState.mt(0) - NoOfajForbk*(returnIfExists(oldfullState.zcoefs, item.a) + returnIfExists(oldfullState.zcoefs, item.b))
        }else{
          SXjkkj - Njkkj * oldfullState.mt(0) - NoOfajForbk*(returnIfExists(oldfullState.zcoefs, item.a) + returnIfExists(oldfullState.zcoefs, item.b)) - NoOfakForbj*(returnIfExists(oldfullState.zcoefs, item.b) + returnIfExists(oldfullState.zcoefs, item.a))
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
        curIndicsEstim(item.a,item.b) = 0.0
        curIndicsEstim(item.b, item.a) = curIndicsEstim(item.a, item.b)
        curThetaEstim(item.a,item.b) = breeze.stats.distributions.Gaussian(info.thetaPriorMean, sqrt(1 / oldfullState.tauabth(1))).draw() // sample from the prior of interactions
        curThetaEstim(item.b, item.a) = curThetaEstim(item.a, item.b)
      }
    })

    oldfullState.copy(thcoefs = curThetaEstim, indics = curIndicsEstim, finalCoefs = curThetaEstim*:*curIndicsEstim)
  }

  override protected def getFileNameToSaveResults(param: String): String = {
    val filePath = getMainFilePath.concat("/try")
    val pathToFiles = Map("mutau" -> filePath.concat("mutau.csv") ,
      "taus" -> filePath.concat("taus.csv"),
      "zetas" -> filePath.concat("zetas.csv"),
      "thetas" -> filePath.concat("thetas.csv"),
      "indics" -> filePath.concat("indics.csv"),
      "allcoefs" -> filePath.concat("allCoefs.csv")
    )
    pathToFiles(param)
  }

  override def getInputFilePath(): String = getMainFilePath.concat("/simulInterSymmetricBoth.csv")

}
