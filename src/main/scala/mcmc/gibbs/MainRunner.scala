package mcmc.gibbs

import java.io.File
import breeze.linalg.{DenseVector, csvread, max}
import structure.{DVStructure, DVStructureIndexedMapMemo}
import scala.io.StdIn.readLine

object MainRunner {
  def main(args: Array[String]): Unit = {

    def setDefaultValues() = {
      val filePath = "/home/antonia/ResultsFromCloud/Report/symmetricNov/asymmetricBoth"
      val inputFile = "/simulInterAsymmetricBoth.csv"
      val outputFile = "/try.csv"
      val outputTimeFile = "/try.txt"
      val caseToRun = "AsymmetricBoth"
      val noOfIterations = 100000
      val thin = 10
      val burnIn = 1000
      val logLikFlag = true

      Arguments(noOfIterations, thin, burnIn, logLikFlag, caseToRun, filePath, inputFile, outputFile, outputTimeFile)
    }

    def setValuesFromArguments(defaultArgs: Arguments, argums: Array[String]) = {
      if(argums.size == 0){
        defaultArgs
      }else if(defaultArgs.getClass.getDeclaredFields.size == args.size){
        Arguments(argums(0).toInt, argums(1).toInt, argums(2).toInt, argums(3).toBoolean, argums(4), argums(5), argums(6), argums(7), argums(8))
        // (noOfIters: Int, thin: Int, burnIn: Int, logLikFlag: Boolean, caseToRun: String, pathToFiles: String, inputFile: String, outputFile: String, outputTimeFile: String)
      }else{
        throw new RuntimeException("Wrong number of arguments passed. \n Arguments need to be: noOfIters: Int, thin: Int, burnIn: Int, logLikFlag: Boolean, caseToRun: String, pathToFiles: String, inputFile: String, outputFile: String, outputTimeFile: String")
      }
    }

    val defaultArgums = setDefaultValues()
    val updatedArgums = setValuesFromArguments(defaultArgums, args)

    val varSelectionObject = getVariableSelectionVariant(updatedArgums.caseToRun)

    varSelectionObject.filesDirectory = updatedArgums.pathToFiles
    varSelectionObject.outputFile = updatedArgums.pathToFiles.concat(updatedArgums.outputFile)
    varSelectionObject.outputTimeFile = updatedArgums.pathToFiles.concat(updatedArgums.outputTimeFile)

    //stop execution until press enter
    //readLine()
    println("Running")

    // Read the data
    val data = csvread(new File(updatedArgums.pathToFiles.concat(updatedArgums.inputFile)))
    val sampleSize = data.rows
    val y = data(::, 0)
    val sumObs = y.toArray.sum // Sum of the values of all the observations
    val alpha = data(::, 1).map(_.toInt).map(x => x - 1)
    val beta = data(::, 2).map(_.toInt).map(x => x - 1)
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

    val structureSorted : DVStructure = new DVStructureIndexedMapMemo(y, alphaSorted, betaSorted) // Sorted structure to be used for the indices to run through the data but not repeat e.g. only (1,3) and not (3,1)
    val noOfInters = structureSorted.sizeOfStructure()
    val zetaLevels = zetaDistinct.max + 1
    val zetaLevelsDist = zetaDistinct.length

    val sizeofDouble = structure.sizeOfDouble()

    val noOftriangular = zetaLevels * (zetaLevels+1) / 2
    //Used for SymmetricMain
    val alphaLevelsDist = alphaDistinct.length
    val betaLevelsDist = betaDistinct.length

    // Parameters
    val aPrior = 1
    val bPrior = 0.0001
    val alphaPriorMean = 0.0
    val betaPriorMean = 0.0
    val mu0 = 0.0
    val tau0 = 0.0001
    val a = 1
    val b = 0.0001
    val interPriorMean = 0.0 //common mean for all the interaction effects
    val pap = 2
    val pbp = 10

    val initialInfo = InitialInfo(updatedArgums.noOfIters, updatedArgums.thin, updatedArgums.burnIn, sampleSize, sumObs, structure, structureSorted, alphaLevels, betaLevels, zetaLevels, noOfInters, sizeofDouble, alphaLevelsDist, betaLevelsDist,  zetaLevelsDist, noOftriangular,
      alphaPriorMean, betaPriorMean, interPriorMean, mu0, tau0,
      a, b, aPrior, bPrior, pap, pbp, updatedArgums.logLikFlag)

    varSelectionObject.time(
      varSelectionObject.variableSelection(initialInfo)
    )
  }

  private def getVariableSelectionVariant(s: String) : VariableSelection = {

    val objectToRun = if (s.equalsIgnoreCase("AsymmetricBoth")){
      object myAsymmetricBoth extends AsymmetricBoth
      myAsymmetricBoth
    } else if (s.equalsIgnoreCase("SymmetricInteractions")){
      object mySymmetricInters extends SymmetricInters
      mySymmetricInters
    }else if (s.equalsIgnoreCase("SymmetricMain")) {
      object mySymmetricMain extends SymmetricMain
      mySymmetricMain
    }else if (s.equalsIgnoreCase("SymmetricBoth")) {
      object mySymmetricBoth extends SymmetricBoth
      mySymmetricBoth
    }else if (s.equalsIgnoreCase("Saturated")) {
      object mySatAsymmetricBoth extends SaturatedAsymmetricBoth
      mySatAsymmetricBoth
    }else{
      throw new RuntimeException("Wrong model declaration. Should be one of: \"AsymmetricBoth\", \"SymmetricInteractions\", \"SymmetricMain\", \"SymmetricBoth\", \"Saturated\" ")
    }
    objectToRun
  }
}

