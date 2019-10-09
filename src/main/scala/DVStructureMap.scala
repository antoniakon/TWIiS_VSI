import breeze.linalg.{DenseVector, max}

class DVStructureMap(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int]) extends DVStructure {

  private val myStructure = scala.collection.mutable.Map[(Int, Int), DVList]()

  val alphaLevels = alpha.toArray.distinct.length
  val betaLevels = beta.toArray.distinct.length
  val zetaLevels = max(alphaLevels, betaLevels)

  init()

  private def init(): Unit = {
    for (i <- 0 until y.length) {
      myStructure.get(alpha(i), beta(i)) match {
        case None    => myStructure += ((alpha(i), beta(i)) -> new DVList())
        case Some(value) => // do nothing
      }

      myStructure((alpha(i), beta(i))).addItem(y(i))
    }
  }

  /**
    * Calculates the sum of the response y for a given alpha
    */
  override def calcAlphaSum(j: Int): Double = {
    var sum = 0.0
    myStructure.foreach( item =>
      if (item._1._1 == j) {
        sum = sum + item._2.sum
      }
    )
    sum
  }

  /**
    * Calculates the sum of the response y for a given beta
    */
  override def calcBetaSum(k: Int): Double = {
    var sum = 0.0
    myStructure.foreach( item =>
      if (item._1._2 == k) {
        sum = sum + item._2.sum
      }
    )
    sum
  }

  /**
    * Calculates the number of the responses y for a given alpha
    */
  override def calcAlphaLength(j: Int): Double = {
    var length = 0.0
    myStructure.foreach( item =>
      if (item._1._1 == j) {
        length = length + item._2.length
      }
    )
    length
  }

  /**
    * Calculates the number of the responses y for a given beta
    */
  override def calcBetaLength(k: Int): Double = {
    var length = 0.0
    myStructure.foreach( item =>
      if (item._1._2 == k) {
        length = length + item._2.length
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

  override def getAllItemsForGivenA(a: Int): List[DVItem] = {
    myStructure
      .filterKeys(key => key._1 == a)
      .map(item => new DVItem(item._1._1, item._1._2, item._2))
      .toList
  }

  override def getAllItemsForGivenB(b: Int): List[DVItem] = {
    myStructure
      .filterKeys(key => key._2 == b)
      .map(item => new DVItem(item._1._1, item._1._2, item._2))
      .toList
  }

  override def getAllItemsMappedByA(): Map[Int, List[DVItem]] = {
    myStructure.keys.map(k => k._1).map(x => (x, getAllItemsForGivenA(x))).toMap
  }

  override def getAllItemsMappedByB(): Map[Int, List[DVItem]] = {
    myStructure.keys.map(k => k._2).map(x => (x, getAllItemsForGivenB(x))).toMap
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
    getAllItemsMappedByZ()
      .filterKeys(k => k==z).flatMap(elem => elem._2)
  }

  override def getAllItemsMappedByZ(): Map[Int, Map[(Int, Int), DVList]] = {
    var myMap = scala.collection.mutable.Map[Int,  Map[(Int, Int), DVList]]()

    for (i <- 0 until zetaLevels) {
      val selectedItems = myStructure.filterKeys(key => (key._1 == i || key._2 == i) && !(key._1 == i && key._2 == i)).toMap
      myMap += (i -> selectedItems)

    }
    myMap.toMap
  }
}
