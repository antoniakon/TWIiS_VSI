import breeze.linalg.DenseVector

class DVStructureMap(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int]) extends DVStructure {

  private val myStructure = scala.collection.mutable.Map[(Int, Int), DVList]()
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
}
