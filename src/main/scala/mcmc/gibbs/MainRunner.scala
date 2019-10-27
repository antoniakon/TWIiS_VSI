package mcmc.gibbs

import java.io.File

import breeze.linalg.csvread
import structure.{DVStructure, DVStructureMap}

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
    val structure : DVStructure = new DVStructureMap(y, alpha, beta)
    val alphaLevels = alpha.toArray.distinct.length
    val betaLevels = beta.toArray.distinct.length

    // Parameters
    val noOfIters = 10
    val thin = 1
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

    val initialInfo = InitialInfo(noOfIters, thin, sampleSize, sumObs, structure, alphaLevels, betaLevels,
      alphaPriorMean, betaPriorMean, interPriorMean, mu0, tau0,
      a, b, aPrior, bPrior, p)

    val statesResults =
      varSelectionObject.time(
        varSelectionObject.variableSelection(initialInfo)
      )
    varSelectionObject.printResults(statesResults)
  }


  private def getVariableSelectionVariant() : VariableSelection = {
    AsymmetricBoth
  }

}

