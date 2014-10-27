package org.hammerlab.guacamole.commands

import org.apache.spark.SparkContext
import org.apache.spark.storage.StorageLevel
import org.hammerlab.guacamole.Common.Arguments.GermlineCallerArgs
import org.hammerlab.guacamole._
import org.hammerlab.guacamole.filters.GenotypeFilter.GenotypeFilterArguments
import org.hammerlab.guacamole.reads.{ MappedRead, Read }
import org.hammerlab.guacamole.variants.Breakpoint
import org.hammerlab.guacamole.windowing.SlidingWindow
import org.kohsuke.args4j.{ Option => Opt }

object GermlineSV {

  protected class Arguments extends GermlineCallerArgs with GenotypeFilterArguments {
    @Opt(name = "--threshold", metaVar = "X", usage = "Make a call if at least X% of reads support it. Default: 20")
    var threshold: Int = 20
  }

  object Caller extends SparkCommand[Arguments] {

    override val name = "germline-sv"
    override val description = "call structural variants by examining paired end reads"

    override def run(args: Arguments, sc: SparkContext): Unit = {

      val filters = Read.InputFilters(mapped = true, nonDuplicate = true)
      val readSet = Common.loadReadsFromArguments(args, sc, filters = filters)
      val mappedReads = readSet.mappedReads

      Common.progress("Loaded %,d mapped non-duplicate MdTag-containing reads into %,d partitions.".format(
        mappedReads.count, mappedReads.partitions.length))

      val loci = Common.loci(args, readSet)
      val lociPartitions = DistributedUtil.partitionLociAccordingToArgs(
        args,
        loci,
        mappedReads
      )

      val threshold = args.threshold / 100.0

      val genotypes = DistributedUtil.windowTaskFlatMapMultipleRDDs(
        Seq(mappedReads),
        lociPartitions,
        0L,
        (task, taskLoci, taskRegionsSeq: Seq[Iterator[MappedRead]]) => {
          DistributedUtil.collectByContig[MappedRead, Breakpoint](
            taskRegionsSeq,
            taskLoci,
            halfWindowSize = 0L,
            (loci, windows) => {
              val lociIterator = loci.iterator
              var lastBreakpoint: Option[Breakpoint] = None
              val builder = Vector.newBuilder[Breakpoint]
              while (SlidingWindow.advanceMultipleWindows(windows, lociIterator, true).isDefined) {
                val window = windows(0)
                // Find a new breakpoint or extend the last one
                val (currentBreakpoint, variants) = discoverBreakpoints(lastBreakpoint,
                  window.currentRegions(),
                  threshold)
                lastBreakpoint = currentBreakpoint
                builder ++= variants
              }
              builder ++= lastBreakpoint
              builder.result.iterator
            }
          )
        }
      )

      genotypes.persist(StorageLevel.MEMORY_ONLY_SER)
      val minReadDepth = args.minReadDepth
      val minAltReadDepth = args.minAlternateReadDepth
      val filteredGenotypes = genotypes.filter(g => g.support < minReadDepth || g.altSupport < minAltReadDepth)
      Common.progress("Computed %,d structural variants".format(genotypes.count))
      filteredGenotypes.map(_.toString).saveAsTextFile(args.variantOutput)

      DelayedMessages.default.print()

    }

    /**
     *
     * Discover a new breakpoint (inversion, translocation or tandem duplication at this locus
     * If a breakpoint was found early that overlaps this postion it is extended
     * If a breakpoint was found early that no longer overlaps this position a variant is produced
     *
     * @param lastBreakpoint current breakpoint that can be extends
     * @param readsAtLocus reads that overlap this locus
     * @param threshold mininum percent of reads needed to support a breakpoint
     * @return
     */
    def discoverBreakpoints(lastBreakpoint: Option[Breakpoint],
                            readsAtLocus: Seq[MappedRead],
                            threshold: Double = 0.3): (Option[Breakpoint], Seq[Breakpoint]) = {

      val readsSupportingBreakpoint = readsAtLocus.filter(
        read => read.inDuplicatedRegion || read.inInvertedRegion || read.inTranslocatedRegion)
      val breakpointRatio = readsSupportingBreakpoint.size.toFloat / readsAtLocus.size

      // Check for minimum threshold of read that support the breakpoint
      if (breakpointRatio > threshold) {
        val newBreakpoint = Breakpoint(readsSupportingBreakpoint)
        lastBreakpoint match {
          case Some(breakpoint) => {
            if (newBreakpoint.overlaps(breakpoint)) {
              // Merge the two breakpoints
              (Some(newBreakpoint.merge(breakpoint)), Seq.empty)
            } else {
              // Emit a new breakpoint and a variant for the last breakpoint
              (Some(newBreakpoint), Seq(breakpoint))
            }
          }
          case None => (Some(newBreakpoint), Seq.empty)
        }
      } else {
        (lastBreakpoint, Seq.empty) // If there is no new breakpoint, keep the last one
      }

    }

  }

}

