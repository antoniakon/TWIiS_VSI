import scala.collection.mutable.ListBuffer

/**
  * Created by Antonia Kontaratou.
  * A DVList list contains observations of the response y.
  */
class DVList {
  val list= new ListBuffer[Double]

  def addItem(item: Double)= {
    list += item
  }

  def sum: Double= {
    list.sum
  }

  def length: Int= {
    list.length
  }

  override def toString = {
    sum.toString
  }

}
