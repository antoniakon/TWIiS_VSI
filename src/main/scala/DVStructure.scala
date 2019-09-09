import breeze.linalg.DenseVector

/**
  * Created by Antonia Kontaratou.
  * Creates the DVStructure of the form: Array[Array[DVList]]
  */
class DVStructure(y: DenseVector[Double], alpha: DenseVector[Int], beta: DenseVector[Int])  {

  val nj = alpha.toArray.distinct.length
  val nk = beta.toArray.distinct.length
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
  def calcAlphaSum(j: Int): Double = {
    var sum = 0.0
    for(k <- 0 until nk){
      sum = sum + myStructure(j)(k).sum
    }
    sum
  }

  /**
    * Calculates the sum of the response y for a given beta
    */
  def calcBetaSum(k: Int): Double = {
    var sum = 0.0
    for(j <- 0 until nj){
      sum = sum + myStructure(j)(k).sum
    }
    sum
  }

  /**
    * Calculates the number of the responses y for a given alpha
    */
  def calcAlphaLength(j: Int): Double = {
    var length = 0.0
    for(k <- 0 until nk){
      length = length + myStructure(j)(k).length
    }
    length
  }

  /**
    * Calculates the number of the responses y for a given beta
    */
  def calcBetaLength(k: Int): Double = {
    var length = 0.0
    for(j <- 0 until nj){
      length = length + myStructure(j)(k).length
    }
    length
  }
  /**
    * Calculates the number of the responses y for a given alpha and a given beta
    */
  def calcUniNatLength(j: Int, k:Int): Double = {
    var length = 0.0
    for(k <- 0 until nk){
      length = length + myStructure(j)(k).length
    }
    length
  }

  /**
    * Returns the DVList for (alpha,beta)
    */
  def getDVList(j: Int, k: Int): DVList = {
    myStructure(j)(k)
  }

}
