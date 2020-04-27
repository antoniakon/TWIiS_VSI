package mcmc.gibbs

import java.io.File
import breeze.linalg.{DenseVector, csvread, max}
import structure.{DVStructure, DVStructureIndexedMap, DVStructureIndexedMapMemo}
import scala.io.StdIn.readLine

object MainRunner {
  def main(args: Array[String]): Unit = {
    val varSelectionObject = getVariableSelectionVariant()

    //stop execution until press enter
    readLine()

    // Read the data
    val data = csvread(new File(varSelectionObject.getInputFilePath()))
    val sampleSize = data.rows
    val y = data(::, 0)
    val sumObs = y.toArray.sum // Sum of the values of all the observations
    val alpha = data(::, 1).map(_.toInt).map(x => x - 1)
    val beta = data(::, 2).map(_.toInt).map(x => x - 1)
    //    val structure : DVStructure = new DVStructureArrays(y, alpha, beta)
    val structure : DVStructure = new DVStructureIndexedMapMemo(y, alpha, beta)
    val alphaDistinct = alpha.toArray.distinct
    val betaDistinct = beta.toArray.distinct
    val zetaDistinct = alphaDistinct.union(betaDistinct).distinct
    val alphaLevels = alphaDistinct.max + 1
    val betaLevels = betaDistinct.max + 1

    // For the symmetric interactions case sort the data, smaller first
    val alphaSorted = DenseVector.zeros[Int](sampleSize)
    val betaSorted = DenseVector.zeros[Int](sampleSize)
    for (i <- 0 until sampleSize) {
      if (alpha(i) <= beta(i)) {
        alphaSorted(i) = alpha(i)
        betaSorted(i) = beta(i)
      } else {
        alphaSorted(i) = beta(i)
        betaSorted(i) = alpha(i)
      }
    }

    val structureSorted : DVStructure = new DVStructureIndexedMap(y, alphaSorted, betaSorted) // Sorted structure to be used for the indices to run through the data but not repeat e.g. only (1,3) and not (3,1)
    val noOfInters = structureSorted.sizeOfStructure()
    val zetaLevels = zetaDistinct.max + 1
    val zetaLevelsDist = zetaDistinct.length

    val sizeofDouble = structure.sizeOfDouble()

    val noOftriangular = zetaLevels * (zetaLevels+1) / 2
    //Used for SymmetricMain
    val alphaLevelsDist = alphaDistinct.length
    val betaLevelsDist = betaDistinct.length

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
    val burnIn = 10

    val initialInfo = InitialInfo(noOfIters, thin, burnIn, sampleSize, sumObs, structure, structureSorted, alphaLevels, betaLevels, zetaLevels, noOfInters, sizeofDouble, alphaLevelsDist, betaLevelsDist,  zetaLevelsDist, noOftriangular,
      alphaPriorMean, betaPriorMean, interPriorMean, mu0, tau0,
      a, b, aPrior, bPrior, p)

    varSelectionObject.time(
      varSelectionObject.variableSelection(initialInfo)
    )
  }

  private def getVariableSelectionVariant() : VariableSelection = {
    object myAsymmetricBoth extends AsymmetricBoth
    object mySymmetricInters extends SymmetricInters
    object mySymmetricMain extends SymmetricMain
    object mySymmetricBoth extends SymmetricBoth
    object mySatAsymmetricBoth extends SaturatedAsymmetricBoth
    mySatAsymmetricBoth
  }
}

