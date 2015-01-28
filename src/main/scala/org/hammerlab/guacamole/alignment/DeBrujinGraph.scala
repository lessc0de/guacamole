package org.hammerlab.guacamole.alignment

import org.hammerlab.guacamole.Bases
import scala.collection.mutable

class DeBrujinGraph(val kmerCounts: mutable.Map[Kmer, Int], kmerSize: Int) {

  type Path = List[Kmer]

  def pruneKmers(minSupport: Int) = {
    kmerCounts
      .filter(_._2 < minSupport)
      .foreach( {case (kmer, count) => kmerCounts.remove(kmer) })
  }

  def allPaths(seed: Kmer,
               minPathLength: Int = 100,
               maxPathLength: Int = 1000,
               maxPaths: Int = 10): List[(Path, Int)] = {
    var paths: List[(Path, Int)] = List.empty[(Path, Int)]
    var frontier: mutable.Stack[Kmer] = mutable.Stack(seed)
    var currentPath: Path = List.empty
    var pathScore: Int = 0
    // explore branches until we accumulate the maximum number of appropriate length paths
    while (frontier.nonEmpty && paths.size < maxPaths) {
      val next = frontier.pop()
      pathScore += kmerCounts(next)
      currentPath = next :: currentPath
      val nextChildren = children(next)
      if (nextChildren.nonEmpty && currentPath.size < maxPathLength) {
        frontier ++= nextChildren
      } else {
        // avoid short branches and repetitive branches
        if (currentPath.size + 1 >= minPathLength && currentPath.size < maxPathLength)
          paths = (currentPath, pathScore) :: paths
        currentPath = List.empty
        pathScore = 0
      }
    }
    paths
  }

  def children(node: Kmer): Seq[Kmer] = {
    node.possibleNext.filter(kmerCounts.contains).toSeq
  }

}

object DeBrujinGraph {
  type Sequence = Seq[Byte]
  def apply(sequences: Seq[Sequence], kmerSize: Int): DeBrujinGraph = {
    val kmerCounts: mutable.Map[Kmer, Int] = mutable.Map.empty[Kmer, Int]

    sequences.filter(Bases.allStandardBases(_))
      .foreach(
        _.sliding(kmerSize)
          .foreach(seq => {
          val kmer = Kmer(seq)
          val count = kmerCounts.getOrElse(kmer, 0)
          kmerCounts.update(kmer, count + 1)
        })
      )

    new DeBrujinGraph(kmerCounts, kmerSize)
  }

  def apply(sequences: Seq[String], kmerSize: Int, isString: Boolean = true): DeBrujinGraph = {
    DeBrujinGraph(sequences.view.map(Bases.stringToBases(_)), kmerSize)
  }

}
