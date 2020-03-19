package structure

import breeze.linalg.{DenseVector, max}

import scala.collection.mutable.ListBuffer

class DVStructureIndexedMapMemo(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int]) extends DVStructure {

  val alphaLevels = alpha.toArray.distinct.length
  val betaLevels = beta.toArray.distinct.length
  val zetaLevels = max(alphaLevels, betaLevels)

  private val myStructure = scala.collection.mutable.Map[(Int, Int), DVList]()
  private var alphaIndices = scala.collection.mutable.Map[Int, ListBuffer[(Int, Int)]]()
  private var betaIndices = scala.collection.mutable.Map[Int, ListBuffer[(Int, Int)]]()
  private var zetaIndices = scala.collection.mutable.Map[Int, ListBuffer[(Int, Int)]]()
  private var zetaIndicesWithoutDoubles = scala.collection.mutable.Map[Int, ListBuffer[(Int, Int)]]()
  private var zetaIndicesDoubles = scala.collection.mutable.Map[Int, ListBuffer[(Int, Int)]]()

  //private val myStructure: TreeMap[(Int, Int), DVList] = initMap()
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

      alphaIndices = alphaIndices.map{case (k,v) => (k, v.distinct)}
      betaIndices = betaIndices.map{case (k,v) => (k, v.distinct)}
      zetaIndices = zetaIndices.map{case (k,v) => (k, v.distinct)}
      zetaIndicesWithoutDoubles = zetaIndices.map{case (k,v) => (k, v.filter(a => a._1!=a._2))}
      zetaIndicesDoubles = zetaIndices.map{case (k,v) => (k, v.filter(a => a._1==a._2))}.filter(v1 => v1._2.nonEmpty) //Includes only the double z without the zs that do not have doubles
      myStructure((curAlpha, curBeta)).addItem(y(i))

    }
  }


  /**
    * Calculates the sum of the response y for a given alpha
    */
  override def calcAlphaSum(j: Int): Double = {
    val sum = alphaIndices(j).map(tuple => myStructure(tuple).sum).sum
    sum
  }


  /**
    * Calculates the sum of the response y for a given beta
    */
  override def calcBetaSum(k: Int): Double = {
    val sum = betaIndices(k).map(tuple => myStructure(tuple).sum).sum
    sum
  }

  /**
    * Calculates the number of the responses y for a given alpha
    */
  override def calcAlphaLength(j: Int): Double = {
    val length = alphaIndices(j).map(tuple => myStructure(tuple).length).sum
    length
  }

  /**
    * Calculates the number of the responses y for a given beta
    */
  override def calcBetaLength(k: Int): Double = {
    val length = betaIndices(k).map(tuple => myStructure(tuple).length).sum
    length
  }

  /**
    * Calculates the number of the responses y for a given alpha and beta
    */
  override def calcAlphaBetaLength(j: Int, k: Int): Double = {
    // Uses Option because if the key (j,k) is not found it throws a java.util.NoSuchElementException
    val lengthMaybe = myStructure.get(j,k)
    val length = lengthMaybe match {
      case Some(dvlist) =>
        dvlist.length
      case None =>
        0
    }
    length
  }

  /**
    * Calculates the sum of the responses y for a given alpha and beta
    */
  override def calcAlphaBetaSum(j: Int, k: Int): Double = {
    // Uses Option because if the key (j,k) is not found it throws a java.util.NoSuchElementException
    val lengthMaybe = myStructure.get(j,k)
    val sum = lengthMaybe match {
      case Some(dvlist) =>
        dvlist.sum
      case None =>
        0
    }
    sum
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

  /**
    * For a given a it returns a list of items  DVItem(a: Int,b: Int,list: DVList)
    */
  override def getAllItemsForGivenA(a: Int): List[DVItem] = {
    val ifaexists = alphaIndices.get(a)
    val res = ifaexists match {
      case Some(tups) =>
        tups.map(item => new DVItem(item._1, item._2, myStructure(item._1, item._2)))
          .toList
      case None => List[DVItem]()
    }
    res
  }

  /**
    * For a given b it returns a list of items  DVItem(a: Int,b: Int,list: DVList)
    */
  override def getAllItemsForGivenB(b: Int): List[DVItem] = {
    val ifbexists = betaIndices.get(b)
    val res = ifbexists match {
      case Some(tups) =>
        tups.map(item => new DVItem(item._1, item._2, myStructure(item._1, item._2)))
          .toList
      case None => List[DVItem]()
    }
    res
  }

  /**
    * It returns a Map[Int, List[DVItem]], where the keys are the a's
    */
  override def getAllItemsMappedByA(): Map[Int, List[DVItem]] = {
    alphaIndices.map(x => (x._1, getAllItemsForGivenA(x._1))).toMap
  }

  /**
    * It returns a Map[Int, List[DVItem]], where the keys are the b's
    */
  override def getAllItemsMappedByB(): Map[Int, List[DVItem]] = {
    betaIndices.map(x => (x._1, getAllItemsForGivenB(x._1))).toMap
  }

  /**
    * Calculates the sum of the response y for a given zeta
    */
  override def calcZetaSum(zj: Int): Double = {
    val sum = zetaIndicesWithoutDoubles(zj).map(tuple => myStructure(tuple).sum).sum
    sum
  }

  /**
    * Calculates the number of the responses y for a given zeta
    */
  override def calcZetaLength(zj: Int): Double = {
    val length = zetaIndicesWithoutDoubles(zj).map(tuple => myStructure(tuple).length).sum
    length
  }

  /**
    * Calculates the sum of the response y for a given zeta
    */
  override def calcDoubleZetaSum(zj: Int): Double = {
    // val sum = zetaIndicesDoubles(zj).map(tuple => myStructure(tuple).sum).sum
    val lengthMaybe = zetaIndicesDoubles.get(zj)
    val sum = lengthMaybe match {
      case Some(listbuf) =>
        listbuf.map(tuple => myStructure(tuple).sum).sum
      case None =>
        0
    }
    sum
  }

  /**
    * Calculates the number of the responses y for a given zeta
    */
  override def calcDoubleZetaLength(zj: Int): Double = {
    //val length = zetaIndicesDoubles(zj).map(tuple => myStructure(tuple).length).sum
    val lengthMaybe = zetaIndicesDoubles.get(zj)
    val length = lengthMaybe match {
      case Some(listbuf) =>
        listbuf.map(tuple => myStructure(tuple).length).sum
      case None =>
        0
    }
    length
  }

  /**
    * Returns a Map[(Int,Int),DVList] with all the cases where zeta is either on the first side or the second without being in both
    */
  override def getAllOtherZetasItemsForGivenZ(z: Int): Map[(Int,Int),DVList] = zetaIndicesWithoutDoubles(z).map(tuple => (tuple, myStructure(tuple))).toMap

  override def getAllDoubleZetasItemsForGivenZ(z: Int): Map[(Int,Int),DVList] = {
    //  zetaIndicesDoubles(z).map(tuple => (tuple, checkIfExists(myStructure.get(tuple)))).toMap
    //zetaIndicesDoubles(z).map(tuple => (tuple, myStructure(tuple))).toMap

    val maybeList = zetaIndicesDoubles.get(z)
    val returnedMap = maybeList match {
      case Some(listbuf) =>
        listbuf.map(tuple => (tuple, myStructure(tuple))).toMap
      case None =>
        Map((0,0) -> new DVList())
    }
    returnedMap
  }


  /**
    * Returns a Map[(Int,Int),DVList] with all the cases where zeta is either on the first side or the second (both sides included)
    */
  override def getZetasItemsForGivenZ(z: Int): Map[(Int,Int),DVList] = zetaIndices(z).map(tuple => (tuple, myStructure(tuple))).toMap

  override def sizeOfStructure():Int = myStructure.keys.size

  override def sizeOfDouble():Int = myStructure.keys.filter(k => (k._1==k._2)).size
}
