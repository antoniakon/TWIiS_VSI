package structure

import java.util.stream.Collectors

import breeze.linalg.{DenseVector, max}

import scala.collection.JavaConverters._
import scala.collection.mutable.ListBuffer

class DVStructureJavaMap(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int]) extends DVStructure {



  val alphaLevels = alpha.toArray.distinct.length
  val betaLevels = beta.toArray.distinct.length
  val zetaLevels = max(alphaLevels, betaLevels)
  private val newStructure = scala.collection.mutable.Map[Int,  Map[(Int, Int), DVList]]()

  private val myStructure: java.util.HashMap[(Int, Int), DVList] = initMap()
  getAllItemsMappedByZ()
  println(newStructure)

  private def initMap(): java.util.HashMap[(Int, Int), DVList] = {
    val tempMap = new java.util.HashMap[(Int, Int), DVList]()

    for (i <- 0 until y.length) {
      val list = new DVList()
      val list1 = tempMap.putIfAbsent((alpha(i), beta(i)), list)
      val listToAdd = if (list1 == null) {
        list
      } else {
        list1
      }

      listToAdd.addItem(y(i))
    }

    tempMap
  }

  /**
    * Calculates the sum of the response y for a given alpha
    */
  override def calcAlphaSum(j: Int): Double = {
    var sum = 0.0
    myStructure.forEach( (key, value) =>
      if (key._1 == j) {
        sum = sum + value.sum
      }
    )
    sum
  }

  /**
    * Calculates the sum of the response y for a given beta
    */
  override def calcBetaSum(k: Int): Double = {
    var sum = 0.0
    myStructure.forEach( (key, value) =>
      if (key._2 == k) {
        sum = sum + value.sum
      }
    )
    sum
  }

  /**
    * Calculates the number of the responses y for a given alpha
    */
  override def calcAlphaLength(j: Int): Double = {
    var length = 0.0
    myStructure.forEach( (key, value) =>
      if (key._1 == j) {
        length = length + value.length
      }
    )
    length
  }

  /**
    * Calculates the number of the responses y for a given beta
    */
  override def calcBetaLength(k: Int): Double = {
    var length = 0.0
    myStructure.forEach( (key, value) =>
      if (key._2 == k) {
        length = length + value.length
      }
    )
    length
  }

  /**
    * Calculates the number of the responses y for a given alpha and beta
    */
  override def calcAlphaBetaLength(j: Int, k: Int): Double = {
    // Uses Option because if the key (j,k) is not found it throws a java.util.NoSuchElementException
    val lengthMaybe = myStructure.get(j,k)

    if (lengthMaybe != null) {
      lengthMaybe.length
    } else {
      0
    }
  }

  /**
    * Returns the DVList for (alpha,beta)
    */
  override def getDVList(j: Int, k: Int): DVList = {
    myStructure.get(j,k)
  }

  override def foreach[U](f: DVItem => U): Unit = {
    myStructure.forEach( (key, value) => f(new DVItem(key._1, key._2, value)) )
  }

  override def getAllItemsForGivenA(a: Int): List[DVItem] = {
//    myStructure.entrySet().stream()
//      .filter(entry => entry.getKey._1 == a)
//      .iterator().asScala
//      .map(entry => new DVItem(entry.getKey._1, entry.getKey._2, entry.getValue))
//      .toList

    val scalaList: ListBuffer[DVItem] = new ListBuffer[DVItem]

    myStructure.entrySet().stream()
      .filter(entry => entry.getKey._1 == a)
      .forEach(entry => scalaList += new DVItem(entry.getKey._1, entry.getKey._2, entry.getValue))

    scalaList.toList
  }

  override def getAllItemsForGivenB(b: Int): List[DVItem] = {
    val scalaList: ListBuffer[DVItem] = new ListBuffer[DVItem]

    val value = myStructure.entrySet().stream()
      .filter(entry => entry.getKey._2 == b)
          .collect(Collectors.toList())
//      .collect(Collectors.toMap(x => x.getKey, y => y.getValue))
//      .collect(Collectors.toList())

    asScalaBuffer(value).map{entry=> new DVItem(entry.getKey._1, entry.getKey._2, entry.getValue)}.toList


    //value.stream().map(entry => new DVItem(1, 1, entry.getValue))
//    val skata = mapAsScalaMap(value)


//    value
//      .forEach(entry => scalaList += new DVItem(entry.getKey._1, entry.getKey._2, entry.getValue))

//    scalaList.toList
  }

  override def getAllItemsMappedByA(): Map[Int, List[DVItem]] = {
    val scalaList: ListBuffer[(Int, List[DVItem])] = new ListBuffer[(Int, List[DVItem])]

    myStructure.keySet().stream()
      .map[Int](k => k._1) //.map(x => (x, getAllItemsForGivenA(x)))
      .forEach(x => scalaList += new Pair(x, getAllItemsForGivenA(x)))

    scalaList.toMap
  }

  override def getAllItemsMappedByB(): Map[Int, List[DVItem]] = {
    val scalaList: ListBuffer[(Int, List[DVItem])] = new ListBuffer[(Int, List[DVItem])]

    myStructure.keySet().stream()
      .map[Int](k => k._2) //.map(x => (x, getAllItemsForGivenB(x))).toMap
      .forEach(x => scalaList += new Pair(x, getAllItemsForGivenB(x)))

    scalaList.toMap
  }

  /**
    * Calculates the sum of the response y for a given zeta, not include the cases where k==j
    */
  override def calcZetaSum(zj: Int): Double = {
    //toDo: check flatten after map instead of flatMap. if it give
    getAllOtherZetasItemsForGivenZ(zj)
      .map(elem => elem._2.sum)
      .reduce(_+_)

  }

  /**
    * Calculates the number of the responses y for a given zeta
    */
  override def calcZetaLength(zj: Int): Double = {
    getAllOtherZetasItemsForGivenZ(zj)
      .map(elem => elem._2.length)
      .reduce(_+_)
  }

  override def getAllOtherZetasItemsForGivenZ(z: Int): Map[(Int,Int),DVList] = {
    newStructure
      .filterKeys(k => k==z).flatMap(elem => elem._2).toMap
  }

  private def getAllItemsMappedByZ(): Unit = {

//    for (i <- 0 until zetaLevels) {
//      val selectedItems = myStructure.filterKeys(key => (key._1 == i || key._2 == i) && !(key._1 == i && key._2 == i)).toMap
//      newStructure.get(i) match {
//        case None => newStructure += (i -> selectedItems)
//        case Some(value) => newStructure += (i -> selectedItems)
//      }
//    }

  }

  override def sizeOfStructure():Int = ??? //myStructure.keys.size

  override def getZetasItemsForGivenZ(z: Int): Map[(Int,Int),DVList] = ???

  override def sizeOfDouble():Int = ??? //myStructure.keys.filter(k => (k._1==k._2)).size
  override def calcDoubleZetaSum(zj: Int): Double = ???

  override def calcDoubleZetaLength(zj: Int): Double = ???

  override def getAllDoubleZetasItemsForGivenZ(z: Int): Map[(Int,Int),DVList] = ???

  override def getAllZetas(): List[Int] = ???
}
