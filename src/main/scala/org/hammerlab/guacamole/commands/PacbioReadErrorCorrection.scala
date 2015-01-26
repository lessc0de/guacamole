package org.hammerlab.guacamole.commands

import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.cli.Args4jBase
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.formats.avro.NucleotideContigFragment
import org.hammerlab.guacamole.{DistributedUtil, SparkCommand}
import org.kohsuke.args4j.{Option => Opt}

object PacbioReadErrorCorrection {

  protected class Arguments extends DistributedUtil.Arguments {
    @Opt(name = "--out")
    var output: String = ""

    @Opt(name = "--kmer-size")
    var kmerSize: Int = 14


    @Opt(name = "--long-reads", metaVar = "X", required = true, usage = "Long reads i.e. PacBio")
    var longReads: String = ""

    @Opt(name = "--short-reads", metaVar = "X", required = true, usage = "Short reads i.e. Illumina")
    var shortReads: String = ""
  }

  object Caller extends SparkCommand[Arguments] {
    override val name = "pb-correct"
    override val description = "Correct errors in PacBio long reads using short reads"

    def kmerHashInSequence(sequence: String, kmerSize: Int): Set[Int] = {
      sequence.sliding(kmerSize).map(_.hashCode()).toSet
    }

    override def run(args: Arguments, sc: SparkContext): Unit = {

      val longReads: RDD[NucleotideContigFragment] = sc.adamLoad(args.longReads)
      val shortReads: RDD[NucleotideContigFragment] = sc.adamLoad(args.shortReads)

      longReads.repartition(args.parallelism)
      // Discover all kmers in the long reads
      // Make a map of kmer to partition(s)
      val kmerPartitions = longReads.mapPartitionsWithIndex((partition: Int, reads: Iterator[NucleotideContigFragment]) => {

        reads.map(read => kmerHashInSequence(read.getFragmentSequence, args.kmerSize)).reduce(_ union _).toIterator.map(kc => (kc, partition))
      }).collect().toMap

      // map short reads to the partitions that have that kmer
      val mappedShortReads: RDD[(NucleotideContigFragment, Int)] = shortReads.flatMap(read => {
        kmerHashInSequence(read.getFragmentSequence, args.kmerSize).map( kc => (read, kmerPartitions(kc)))
      })

      // shuffle and partition mappedShortReads
//      val mixedReads = longReads.zipPartitions(mappedShortReads)(
//        (longReadsP: Iterator[NucleotideContigFragment], shortReadsP: Iterator[(NucleotideContigFragment, Int)])  =>
//          longReadsP.foreach( lr => {
//             //map short read to long read )
//          })
//      )
    }
  }
}