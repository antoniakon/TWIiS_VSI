package mcmc.gibbs

import java.io._

import breeze.linalg.{*, DenseMatrix, csvread}
import breeze.stats.mean

abstract class VariableSelection {
  def getInputFilePath(): String
  def variableSelection(info: InitialInfo): FullStateList
  protected def nextmutau(oldfullState: FullState, info: InitialInfo): FullState
  protected def nexttaus(oldfullState: FullState, info: InitialInfo):FullState
  protected def nextCoefs(oldfullState: FullState, info: InitialInfo):FullState
  protected def nextIndicsInters(oldfullState: FullState, info: InitialInfo):FullState
  def printResults(statesResults: FullStateList): Unit
  protected def getFileNameToSaveResults(param: String): String

  @annotation.tailrec
  protected final def calculateNewState(n:Int, info: InitialInfo, fstate:FullState, fstateList:FullStateList): FullStateList = {
    if (n==0) fstateList
    else{
      println(n)
      val latestmt = nextmutau(fstate, info)
      val latesttaus = nexttaus(latestmt, info)
      val latestcoefs = nextCoefs(latesttaus, info)
      val latestFullyUpdatedState = nextIndicsInters(latestcoefs, info)
      if((n % info.thin).equals(0)) {
        calculateNewState(n-1, info, latestFullyUpdatedState, FullStateList(latestFullyUpdatedState::fstateList.fstateL))
      }
      else calculateNewState(n-1, info, latestFullyUpdatedState, fstateList)
    }
  }

  protected final def calculateAndPrintCommons(statesResults: FullStateList): List[DenseMatrix[Double]] = {
    println("mu, tau")
    val mtcoefficients = statesResults.fstateL.map(f => f.mt)
    val mtcoefMat = DenseMatrix(mtcoefficients.map(_.toArray): _*)
    val meanValsmtcoef = mean(mtcoefMat(::, *))
    println(meanValsmtcoef)

    println("thetas")
    val thetascoefficients = statesResults.fstateL.map(f => f.thcoefs.toDenseVector)
    val thetascoefMat = DenseMatrix(thetascoefficients.map(_.toArray): _*)
    val meanValsthetascoef = mean(thetascoefMat(::, *))
    println(meanValsthetascoef)

    println("indicators")
    val indicscoefficients = statesResults.fstateL.map(f => f.indics.toDenseVector)
    val indicscoefMat = DenseMatrix(indicscoefficients.map(_.toArray): _*)
    val meanValsindicscoef = mean(indicscoefMat(::, *))
    println(meanValsindicscoef)

    println("finalCoefs")
    val finalcoefficients = statesResults.fstateL.map(f => f.finalCoefs.toDenseVector)
    val finalcoefMat = DenseMatrix(finalcoefficients.map(_.toArray): _*)
    val meanValsfinalcoef = mean(finalcoefMat(::, *))
    println(meanValsfinalcoef)

    println("taus")
    val tauscoefficients = statesResults.fstateL.map(f => f.tauabth)
    val tauscoefMat = DenseMatrix(tauscoefficients.map(_.toArray): _*)
    val meanValstauscoef = mean(tauscoefMat(::, *))
    println(meanValstauscoef)

    List(mtcoefMat, tauscoefMat, finalcoefMat, indicscoefMat)
  }

  protected final def saveToCSV(mergedMatrix: DenseMatrix[Double], filename: String): Unit = {
    val outputFile = new File(filename)
    breeze.linalg.csvwrite(outputFile, mergedMatrix, separator = ',')
  }

  protected def getMainFilePath() : String ={
    "/home/antonia/ResultsFromCloud/Report/symmetricNov/asymmetricBoth"
  }

  // Calculation of the execution time
  final def time[A](f: => A): A = {
    val s = System.nanoTime
    val ret = f
    val execTime = (System.nanoTime - s) / 1e6
    println("time: " + execTime + "ms")
    val file = new File(getMainFilePath.concat("/ScalaRuntime10mAsymmetricBoth.txt"))
    val bw = new BufferedWriter(new FileWriter(file))
    bw.write(execTime.toString)
    bw.close()
    ret
  }

}
