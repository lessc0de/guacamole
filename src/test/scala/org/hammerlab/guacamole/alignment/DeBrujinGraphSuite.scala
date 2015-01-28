package org.hammerlab.guacamole.alignment

import org.hammerlab.guacamole.TestUtil
import org.hammerlab.guacamole.TestUtil.SparkFunSuite

class DeBrujinGraphSuite extends SparkFunSuite {

  lazy val reads = TestUtil.loadReads(sc, "NA12878_S1-chr1-10000.sam").mappedReads.collect()
  lazy val smallWindowReads = reads.filter(read => read.start > 10000 && read.start < 10100)

  lazy val smallWindowSequences = {
    smallWindowReads.map(_.sequence)
  }

  test("build graph") {

    val sequence = "TCATCTCAAAAGAGATCGA"
    val graph = DeBrujinGraph(Seq(sequence), kmerSize = 8)

    val firstKmer = Kmer("TCATCTCA")
    val nextKmer = Kmer("CATCTCAA")
    val lastKmer = Kmer("GAGATCGA")

    assert(graph.kmerCounts.contains(firstKmer))
    assert(graph.kmerCounts.contains(nextKmer))
    assert(graph.kmerCounts.contains(lastKmer))
  }


  test("build graph with short kmers and correct counts") {

    val sequence = "TCATCTTAAAAGACATAAA"
    val graph = DeBrujinGraph(Seq(sequence), kmerSize = 3)

    val firstKmer = Kmer("TCA")
    val nextKmer = Kmer("CAT")
    val lastKmer = Kmer("AAA")

    assert(graph.kmerCounts(firstKmer) === 1)
    assert(graph.kmerCounts(nextKmer) === 2)
    assert(graph.kmerCounts(lastKmer) === 3)
  }

  test("find all paths") {

    val sequence = "TCATCTTAAAAGACATAAA"
    val graph = DeBrujinGraph(Seq(sequence), kmerSize = 6)
    val paths = graph.allPaths(Kmer("TCATCT"), minPathLength = 10)

    assert(paths.length === 1)
  }

  sparkTest("simple read assembly") {
    val graph = DeBrujinGraph(smallWindowSequences, kmerSize = 20)
    val paths = graph.allPaths(
      Kmer("TAACCCTAACCCTAACCCTA"),
      minPathLength = 85,
      maxPathLength = 125
    )
    assert (paths.length === 5)
  }
}
