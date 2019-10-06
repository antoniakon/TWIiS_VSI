trait DVStructure {

  /**
    * Calculates the sum of the response y for a given alpha
    */
  def calcAlphaSum(j: Int): Double

  /**
    * Calculates the sum of the response y for a given beta
    */
  def calcBetaSum(k: Int): Double

  /**
    * Calculates the number of the responses y for a given alpha
    */
  def calcAlphaLength(j: Int): Double

  /**
    * Calculates the number of the responses y for a given beta
    */
  def calcBetaLength(k: Int): Double

  /**
    * Returns the DVList for (alpha,beta)
    */
  def getDVList(j: Int, k: Int): DVList

  def foreach[U](f: DVItem => U): Unit

  def getAllItemsForGivenA(a : Int): List[DVItem]

  def getAllItemsForGivenB(b : Int): List[DVItem]

  def getAllItemsMappedByA() : Map[Int, List[DVItem]]

  def getAllItemsMappedByB() : Map[Int, List[DVItem]]


}
