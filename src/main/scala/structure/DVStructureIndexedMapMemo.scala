package structure

import breeze.linalg.{DenseVector, max}
import scala.collection.mutable.ListBuffer
import scalaz.Memo

class DVStructureIndexedMapMemo(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int]) extends DVStructure {

  private val EMPTY_DVLIST = new DVList()
  private val EMPTY_LISTBUFF = new ListBuffer[(Int, Int)]()

  val alphaLevels = alpha.toArray.distinct.max+1
  val betaLevels = beta.toArray.distinct.max+1
  val zetaLevels = max(alphaLevels, betaLevels)

  private val myStructure = scala.collection.mutable.Map[(Int, Int), DVList]()
  private var alphaIndices = scala.collection.mutable.Map[Int, ListBuffer[(Int, Int)]]()
  private var betaIndices = scala.collection.mutable.Map[Int, ListBuffer[(Int, Int)]]()
  private var zetaIndices = scala.collection.mutable.Map[Int, ListBuffer[(Int, Int)]]()
  private var zetaIndicesWithoutDoubles = scala.collection.mutable.Map[Int, ListBuffer[(Int, Int)]]()
  private var zetaIndicesDoubles = scala.collection.mutable.Map[Int, ListBuffer[(Int, Int)]]()
  private var allItemsMappedbyA = Map[Int, List[DVItem]]()
  private var allItemsMappedbyB = Map[Int, List[DVItem]]()
  private val alphaBetaLengthMat = Array.ofDim[Double](zetaLevels, zetaLevels).map(x=> x.map(y => 0.0))

  init()

  private def init(): Unit = {
    for (i <- 0 until y.length) {
      val curAlpha = alpha(i)
      val curBeta = beta(i)

      alphaIndices.get(curAlpha) match {
        case None => alphaIndices += curAlpha -> ListBuffer[(Int, Int)]((curAlpha, curBeta))
        case Some(value) => value += ((curAlpha, curBeta)) //could already exist but later use distinct so as not to have duplicates
      }

      betaIndices.get(curBeta) match {
        case None => betaIndices += curBeta -> ListBuffer[(Int, Int)]((curAlpha, curBeta))
        case Some(value) => value += ((curAlpha, curBeta)) //could already exist but later use distinct so as not to have duplicates
      }

      zetaIndices.get(curAlpha) match {
        case None => zetaIndices += curAlpha -> ListBuffer[(Int, Int)]( (curAlpha, curBeta) )
        case Some(value) => value += ( (curAlpha, curBeta) ) //could already exist but later use distinct so as not to have duplicates
      }

      zetaIndices.get(curBeta) match {
        case None => zetaIndices += curBeta -> ListBuffer[(Int, Int)]( (curAlpha, curBeta) )
        case Some(value) => value += ( (curAlpha, curBeta) ) //could already exist but later use distinct so as not to have duplicates
      }

      myStructure.get(curAlpha, curBeta) match {
        case None    => myStructure += ((curAlpha, curBeta) -> new DVList())
        case Some(value) => // do nothing
      }

      myStructure((curAlpha, curBeta)).addItem(y(i))
    }

    alphaIndices = alphaIndices.map{case (k,v) => (k, v.distinct)}
    betaIndices = betaIndices.map{case (k,v) => (k, v.distinct)}
    zetaIndices = zetaIndices.map{case (k,v) => (k, v.distinct)}
    zetaIndicesWithoutDoubles = zetaIndices.map{case (k,v) => (k, v.filter(a => a._1!=a._2))}
    zetaIndicesDoubles = zetaIndices.map{case (k,v) => (k, v.filter(a => a._1==a._2))}.filter(v1 => v1._2.nonEmpty) //Includes only the double z without the zs that do not have doubles

    myStructure.foreach { case (key, value) =>
      alphaBetaLengthMat(key._1)(key._2) = value.length
    }
  }

  private val memoizedCalcAlphaSum: Int => Double = Memo.immutableHashMapMemo {
    num => alphaIndices(num).map(tuple => myStructure(tuple).sum).sum
  }

  /**
    * Calculates the sum of the response y for a given alpha
    */
  override def calcAlphaSum(j: Int): Double = {
    memoizedCalcAlphaSum(j)
  }

  private val memoizedCalcBetaSum: Int => Double = Memo.immutableHashMapMemo {
    num => betaIndices(num).map(tuple => myStructure(tuple).sum).sum
  }

  /**
    * Calculates the sum of the response y for a given beta
    */
  override def calcBetaSum(k: Int): Double = {
    memoizedCalcBetaSum(k)
  }

  private val memoizedCalcAlphaLength: Int => Double = Memo.immutableHashMapMemo {
    num => alphaIndices(num).map(tuple => myStructure(tuple).length).sum
  }
  /**
    * Calculates the number of the responses y for a given alpha
    */
  override def calcAlphaLength(j: Int): Double = {
    memoizedCalcAlphaLength(j)
  }

  private val memoizedCalcBetaLength: Int => Double = Memo.immutableHashMapMemo {
    num => betaIndices(num).map(tuple => myStructure(tuple).length).sum
  }
  /**
    * Calculates the number of the responses y for a given beta
    */
  override def calcBetaLength(k: Int): Double = {
    memoizedCalcBetaLength(k)
  }

  private val memoizedcalcAlphaBetaLength:  Tuple2[Int, Int] => Double = Memo.immutableHashMapMemo {
    // Uses Option because if the key (j,k) is not found it throws a java.util.NoSuchElementException
    num => {
      myStructure.getOrElse((num._1, num._2), EMPTY_DVLIST).length
    }
  }

  /**
    * Calculates the number of the responses y for a given alpha and beta
    */
  override def calcAlphaBetaLength(j: Int, k: Int): Double = {
    //memoizedcalcAlphaBetaLength(j,k)
    alphaBetaLengthMat(j)(k)
  }

  /**
    * Returns the DVList for (alpha,beta)
    */
  override def getDVList(j: Int, k: Int): DVList = {
    myStructure(j,k)
  }

  override def foreach[U](f: DVItem => U): Unit = {
    myStructure.foreach( item => f(new DVItem(item._1._1, item._1._2, item._2)) )
  }

  private val memoizedgetAllItemsForGivenA: Int => List[DVItem] = Memo.immutableHashMapMemo {
    num => {
      alphaIndices.getOrElse(num, EMPTY_LISTBUFF)
        .map(item => new DVItem(item._1, item._2, myStructure(item._1, item._2)))
        .toList
    }
  }
  /**
    * For a given a it returns a list of items  DVItem(a: Int,b: Int,list: DVList)
    */
  override def getAllItemsForGivenA(a: Int): List[DVItem] = {
    memoizedgetAllItemsForGivenA(a)
  }

  private val memoizedgetAllItemsForGivenB: Int => List[DVItem] = Memo.immutableHashMapMemo {
    num => {
      betaIndices.getOrElse(num, EMPTY_LISTBUFF)
        .map(item => new DVItem(item._1, item._2, myStructure(item._1, item._2)))
        .toList
    }
  }

  /**
    * For a given b it returns a list of items  DVItem(a: Int,b: Int,list: DVList)
    */
  override def getAllItemsForGivenB(b: Int): List[DVItem] = {
    memoizedgetAllItemsForGivenB(b)
  }

  allItemsMappedbyA = alphaIndices.map(x => (x._1, getAllItemsForGivenA(x._1))).toMap
  allItemsMappedbyB = betaIndices.map(x => (x._1, getAllItemsForGivenB(x._1))).toMap
  /**
    * It returns a Map[Int, List[DVItem]], where the keys are the a's
    */
  override def getAllItemsMappedByA(): Map[Int, List[DVItem]] = {
    allItemsMappedbyA
  }

  /**
    * It returns a Map[Int, List[DVItem]], where the keys are the b's
    */
  override def getAllItemsMappedByB(): Map[Int, List[DVItem]] = {
    allItemsMappedbyB
  }

  private val memoizedCalcZetaSum: Int => Double = Memo.immutableHashMapMemo {
    num => zetaIndicesWithoutDoubles(num).map(tuple => myStructure(tuple).sum).sum
  }

  /**
    * Calculates the sum of the response y for a given zeta
    */
  override def calcZetaSum(zj: Int): Double = {
    memoizedCalcZetaSum(zj)
  }

  private val memoizedCalcZetaLength: Int => Double = Memo.immutableHashMapMemo {
    num => zetaIndicesWithoutDoubles(num).map(tuple => myStructure(tuple).length).sum
  }

  /**
    * Calculates the number of the responses y for a given zeta
    */
  override def calcZetaLength(zj: Int): Double = {
    memoizedCalcZetaLength(zj)
  }

  private val memoizedCalcDoubleZetaSum: Int => Double = Memo.immutableHashMapMemo {
    num => {
      zetaIndicesDoubles.getOrElse(num, EMPTY_LISTBUFF)
        .map(tuple => myStructure(tuple).sum)
        .sum
    }
  }

  /**
    * Calculates the sum of the response y for a given zeta
    */
  override def calcDoubleZetaSum(zj: Int): Double = {
    memoizedCalcDoubleZetaSum(zj)
  }

  private val memoizedCalcDoubleZetaLength: Int => Double = Memo.immutableHashMapMemo {
    num => {
      zetaIndicesDoubles.getOrElse(num, EMPTY_LISTBUFF)
        .map(tuple => myStructure(tuple).length)
        .sum
    }
  }

  /**
    * Calculates the number of the responses y for a given zeta
    */
  override def calcDoubleZetaLength(zj: Int): Double = {
    memoizedCalcDoubleZetaLength(zj)
  }

  private val memoizedAllOtherZetasItemsForGivenZ: Int => Map[(Int,Int),DVList] = Memo.immutableHashMapMemo {
    num => zetaIndicesWithoutDoubles(num).map(tuple => (tuple, myStructure(tuple))).toMap
  }
  /**
    * Returns a Map[(Int,Int),DVList] with all the cases where zeta is either on the first side or the second without being in both
    */
  override def getAllOtherZetasItemsForGivenZ(z: Int): Map[(Int,Int),DVList] = memoizedAllOtherZetasItemsForGivenZ(z)

  private val memoizedAllDoubleZetasItemsForGivenZ: Int => Map[(Int,Int),DVList] = Memo.immutableHashMapMemo {
    num => {
      zetaIndicesDoubles.getOrElse(num, EMPTY_LISTBUFF)
        .map(tuple => (tuple, myStructure(tuple)))
        .toMap
    }
  }

  override def getAllDoubleZetasItemsForGivenZ(z: Int): Map[(Int,Int),DVList] = {
    memoizedAllDoubleZetasItemsForGivenZ(z)
  }

  private val memoizedZetasItemsForGivenZ: Int => Map[(Int,Int),DVList] = Memo.immutableHashMapMemo {
    num => zetaIndices(num).map(tuple => (tuple, myStructure(tuple))).toMap
  }

  override def getAllZetas(): List[Int] =
    zetaIndices.keys.toList
  /**
    * Returns a Map[(Int,Int),DVList] with all the cases where zeta is either on the first side or the second (both sides included)
    */
  override def getZetasItemsForGivenZ(z: Int): Map[(Int,Int),DVList] = memoizedZetasItemsForGivenZ(z)

  override def sizeOfStructure():Int = myStructure.keys.size

  override def sizeOfDouble():Int = myStructure.keys.filter(k => (k._1==k._2)).size
}
