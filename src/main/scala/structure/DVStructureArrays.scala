package structure

import breeze.linalg.DenseVector

import scala.collection.mutable.ListBuffer

/**
  * Created by Antonia Kontaratou.
  * Creates the DVStructure of the form: Array[Array[DVList]]
  */
class DVStructureArrays(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int]) extends DVStructure {

  val nj = alpha.toArray.distinct.length
  val nk = beta.toArray.distinct.length
  //toDo: change the value of nz
  private val nz = nk
  private val myStructure = Array.ofDim[DVList](nj, nk)
  init

  // Initialise values Array[Array[DVList]]

  /**
    * Initialise values Array[Array[DVList]]
    */
  private def init = {
    for (j <- 0 to nj - 1) {
      for (k <- 0 to nk - 1) {
        myStructure(j)(k) = new DVList()
      }
    }
    for (i <- 0 until y.length) {
      myStructure(alpha(i))(beta(i)).addItem(y(i))
    }
  }

  /**
    * Calculates the sum of the response y for a given alpha
    */
  override def calcAlphaSum(j: Int): Double = {
    var sum = 0.0
    for(k <- 0 until nk){
      sum = sum + myStructure(j)(k).sum
    }
    sum
  }

  /**
    * Calculates the sum of the response y for a given zeta, not include the cases where k==j
    */
  override def calcZetaSum(j: Int): Double = {
    var sum = 0.0
    for(k <- 0 until nk){
      sum = sum + myStructure(j)(k).sum
    }
    sum - myStructure(j)(j).sum
  }

  /**
    * Calculates the sum of the response y for a given beta
    */
  override def calcBetaSum(k: Int): Double = {
    var sum = 0.0
    for(j <- 0 until nj){
      sum = sum + myStructure(j)(k).sum
    }
    sum
  }

  /**
    * Calculates the number of the responses y for a given alpha
    */
  override def calcAlphaLength(j: Int): Double = {
    var length = 0.0
    for(k <- 0 until nk){
      length = length + myStructure(j)(k).length
    }
    length
  }

  /**
    * Calculates the number of the responses y for a given alpha
    */
  override def calcZetaLength(j: Int): Double = {
    var length = 0.0
    for(k <- 0 until nk){
      length = length + myStructure(j)(k).length
    }
    length - myStructure(j)(j).length
  }

  /**
    * Calculates the number of the responses y for a given beta
    */
  override def calcBetaLength(k: Int): Double = {
    var length = 0.0
    for(j <- 0 until nj){
      length = length + myStructure(j)(k).length
    }
    length
  }

  /**
    * Returns the DVList for (alpha,beta)
    */
  override def getDVList(j: Int, k: Int): DVList = {
    myStructure(j)(k)
  }

  override def foreach[U](f: DVItem => U): Unit = {
    for (j <- 0 to nj - 1) {
      for (k <- 0 to nk - 1) {
        f(new DVItem(j, k, myStructure(j)(k)))
      }
    }
  }

  override def getAllItemsForGivenA(a : Int): List[DVItem] = {
    val listBuffer = ListBuffer[DVItem]()

    for ( i <- 0 until nk) {
      listBuffer +=new DVItem(a, i, myStructure(a)(i))
    }

    listBuffer.toList
  }

  override def getAllItemsForGivenB(b : Int): List[DVItem] = {
    val listBuffer = ListBuffer[DVItem]()

    for ( i <- 0 until nj) {
      listBuffer += new DVItem(i, b, myStructure(i)(b))
    }

    listBuffer.toList
  }

  override def getAllItemsMappedByA() : Map[Int, List[DVItem]] = (0 until nj).map(x => (x, getAllItemsForGivenA(x))).toMap

  override def getAllOtherZetasItemsForGivenZ(z: Int): Map[(Int,Int),DVList]= ???

  override def calcAlphaBetaLength(j: Int, k: Int): Double = ???

  override def calcAlphaBetaSum(j: Int, k: Int): Double = ???

  override def sizeOfStructure():Int = ???

  override def sizeOfDouble():Int = ???

  override def getZetasItemsForGivenZ(z: Int): Map[(Int,Int),DVList] = ???

  override def calcDoubleZetaSum(zj: Int): Double = ???

  override def calcDoubleZetaLength(zj: Int): Double = ???

  override def getAllDoubleZetasItemsForGivenZ(z: Int): Map[(Int,Int),DVList] = ???

  override def getAllItemsMappedByB() : Map[Int, List[DVItem]] = (0 until nk).map(x => (x, getAllItemsForGivenB(x))).toMap

  override def getAllZetas(): List[Int] = ???
}
