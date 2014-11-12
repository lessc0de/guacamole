package org.bdgenomics.guacamole.benchmarks

import org.bdgenomics.guacamole.TestUtil
import org.bdgenomics.guacamole.TestUtil.HasSparkContext
import org.bdgenomics.guacamole.pileup.Pileup
import org.scalameter.api._

object PileupBenchmark extends PerformanceTest.Quickbenchmark with HasSparkContext {

  val loci: Gen[Int] = Gen.range("loci")(5, 50, 10)

  val ranges = for {
    size <- loci
  } yield 1 until size


  createSpark("PileupBenchmark", true)

  def loadPileup(filename: String, locus: Long = 0): Pileup = {
    val records = TestUtil.loadReads(sc, filename).mappedReads
    val localReads = records.collect
    Pileup(localReads, locus)
  }

  val pileup = loadPileup("same_start_reads.sam", 0)

  performance of "Pileup" in {
    measure method "advanceToLocus" in {
      using(ranges) in {
        loci =>  loci.foreach(pileup.atGreaterLocus(_, Seq.empty.iterator))
      }
    }
  }
}
