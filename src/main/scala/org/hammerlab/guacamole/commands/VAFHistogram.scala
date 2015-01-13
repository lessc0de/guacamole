package org.hammerlab.guacamole.commands

import breeze.linalg.DenseVector
import com.esotericsoftware.kryo.io.Input
import com.esotericsoftware.kryo.{ Kryo, Serializer, io }
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics
import org.apache.spark.SparkContext
import org.apache.spark.mllib.clustering.{ GaussianMixtureModel, GaussianMixtureEM }
import org.apache.spark.mllib.linalg.Vectors
import org.apache.spark.rdd.RDD
import org.apache.spark.storage.StorageLevel
import org.hammerlab.guacamole.Common.Arguments.{ Output, TumorNormalReads }
import org.hammerlab.guacamole.filters.PileupFilter.PileupFilterArguments
import org.hammerlab.guacamole.pileup.Pileup
import org.hammerlab.guacamole.reads.{ MappedRead, Read }
import org.hammerlab.guacamole.{ Common, DistributedUtil, LociMap, SparkCommand }
import org.kohsuke.args4j.{ Option => Opt }

case class VariantLocus(locus: Long, variantAlleleFrequency: Double)

class VariantLocusSerializer extends Serializer[VariantLocus] {
  override def write(kryo: Kryo, output: io.Output, obj: VariantLocus) = {
    output.writeLong(obj.locus, true)
    output.writeDouble(obj.variantAlleleFrequency)
  }

  override def read(kryo: Kryo, input: Input, clazz: Class[VariantLocus]): VariantLocus = {
    val locus = input.readLong(true)
    val vaf = input.readDouble()

    VariantLocus(locus, vaf)

  }
}

object MixtureCaller {

  protected class Arguments extends DistributedUtil.Arguments with Output with PileupFilterArguments with TumorNormalReads {
    @Opt(name = "--readsToSample", usage = "Minimum log odds threshold for possible variant candidates")
    var readsToSample: Int = 1000
  }

  object Caller extends SparkCommand[Arguments] {
    override val name = "somatic-mixture"
    override val description = "cluster variant allele frequencies"

    override def run(args: Arguments, sc: SparkContext): Unit = {

      val filters = Read.InputFilters(mapped = true, nonDuplicate = true, passedVendorQualityChecks = true)
      val (tumorReads, normalReads) = Common.loadTumorNormalReadsFromArguments(args, sc, filters)

      assert(tumorReads.sequenceDictionary == normalReads.sequenceDictionary,
        "Tumor and normal samples have different sequence dictionaries. Tumor dictionary: %s.\nNormal dictionary: %s."
          .format(tumorReads.sequenceDictionary, normalReads.sequenceDictionary))

      val readsToSample = args.readsToSample

      val loci = Common.loci(args, normalReads)
      val lociPartitions = DistributedUtil.partitionLociAccordingToArgs(
        args,
        loci,
        tumorReads.mappedReads,
        normalReads.mappedReads

      )

      val normalVAFs = generateVariantLociFromReads(normalReads.mappedReads, lociPartitions)
      val tumorVAFs = generateVariantLociFromReads(tumorReads.mappedReads, lociPartitions)

      val tumorVAFHistogram = generateVAFHistogram(tumorVAFs)
      val normalVAFHistogram = generateVAFHistogram(normalVAFs)

      (1 until 99).foreach(
        vaf =>
          println(f"VAF = $vaf - num tumor locus = ${tumorVAFHistogram(vaf)} - num normal locus = ${normalVAFHistogram(vaf)}"))

      buildMixtureModel(tumorVAFs)

      buildMixtureModel(normalVAFs)

    }

    def generateVAFHistogram(vafs: RDD[VariantLocus]): scala.collection.Map[Int, Long] = {
      vafs.map(vaf => math.round(vaf.variantAlleleFrequency * 100).toInt).countByValue()
    }
  }

  def buildMixtureModel(vafs: RDD[VariantLocus]): GaussianMixtureModel = {
    val vafVectors = vafs.map(vaf => Vectors.dense(vaf.variantAlleleFrequency))
    val model = new GaussianMixtureEM()
      .setK(4)
      .setConvergenceTol(1e-2)
      .setMaxIterations(50)
      .run(vafVectors)

    for (i <- 0 until model.k) {
      println(s"Cluster weight=${model.weight(i)}, mean=${model.mu(i)}, spread=${model.sigma(i)}")

    }

    model
  }

  def generateVariantLociFromReads(reads: RDD[MappedRead],
                                   lociPartitions: LociMap[Long],
                                   readsToSample: Int = 1000): RDD[VariantLocus] = {
    val sampleName = reads.take(1)(0).sampleName

    val variantLoci = DistributedUtil.pileupFlatMap[VariantLocus](
      reads,
      lociPartitions,
      skipEmpty = true,
      pileup => generateVariantLocus(pileup).iterator
    )

    variantLoci.persist(StorageLevel.MEMORY_ONLY)
    val numVariantLoci = variantLoci.count
    Common.progress("%d non-zero variant loci in sample %s".format(numVariantLoci, sampleName))

    val sampledVAFs =
      if (numVariantLoci > readsToSample)
        variantLoci
          .sample(withReplacement = false, fraction = readsToSample.toFloat / numVariantLoci)
          .collect()
      else
        variantLoci.collect()

    val stats = new DescriptiveStatistics()
    sampledVAFs.foreach(v => stats.addValue(v.variantAlleleFrequency))

    Common.progress("Variant loci stats (min: %f, max: %f, median: %f, mean: %f, 25: %f, 75: %f)".format(
      stats.getMin, stats.getMax, stats.getPercentile(50), stats.getMean, stats.getPercentile(25), stats.getPercentile(75)
    ))

    variantLoci
  }

  def generateVariantLocus(pileup: Pileup): Option[VariantLocus] = {
    if (pileup.referenceDepth != pileup.depth) {
      Some(VariantLocus(pileup.locus, (pileup.depth - pileup.referenceDepth).toFloat / pileup.depth))
    } else {
      None
    }
  }
}

