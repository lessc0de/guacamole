/**
 * Licensed to Big Data Genomics (BDG) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The BDG licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.hammerlab.guacamole.windowing

import org.apache.spark.Logging
import org.hammerlab.guacamole.DistributedUtil.PerSample
import org.hammerlab.guacamole.{ LociSet, HasReferenceRegion }

import scala.collection.mutable

/**
 * Suppose we have a set of loci on a given contig and some objects that are mapped to that contig, and at each locus we
 * want to look at the regions that overlap a window of a certain number of bases surrounding that locus. This class
 * implements this "sliding window" functionality.
 *
 * After instantiating this class, call [[setCurrentLocus]] repeatedly for each locus being considered. After calling
 * this method, the [[currentRegions]] property will contain the objects that overlap the current locus.
 *
 * To enable an efficient implementation, we require that both the sequence of loci to be considered and the iterator
 * of objects are sorted.
 *
 * @param halfWindowSize Number of nucleotide bases to either side of the specified locus to provide regions for. For
 *                       example, if [[halfWindowSize]]=5, and our [[currentLocus]]=100, then [[currentRegions]] will
 *                       include regions that map to anywhere between 95 and 105, inclusive. Set to 0 to consider only
 *                       regions that overlap the exact locus being considered, with no surrounding window.
 *
 * @param rawSortedRegions Iterator of regions, sorted by the aligned start locus.
 */
case class SlidingWindow[Region <: HasReferenceRegion](halfWindowSize: Long,
                                                       rawSortedRegions: Iterator[Region]) extends Logging {
  /** The locus currently under consideration. */
  var currentLocus = -1L

  /** The new regions that were added to currentRegions as a result of the most recent call to setCurrentLocus. */
  var newRegions: Seq[Region] = Seq.empty

  private var referenceName: Option[String] = None
  private var mostRecentRegionStart: Long = 0
  private val sortedRegions: BufferedIterator[Region] = rawSortedRegions.map(region => {
    if (referenceName.isEmpty) referenceName = Some(region.referenceContig)
    require(region.referenceContig == referenceName.get, "Regions must have the same reference name")
    require(region.start >= mostRecentRegionStart, "Regions must be sorted by start locus")
    mostRecentRegionStart = region.start

    region
  }).buffered

  private val currentRegionsPriorityQueue = {
    // Order regions by end locus, increasing.
    def regionOrdering = new Ordering[Region] {
      def compare(first: Region, second: Region) = second.end.compare(first.end)
    }
    new mutable.PriorityQueue[Region]()(regionOrdering)
  }

  /** The regions that overlap the window surrounding [[currentLocus]]. */
  def currentRegions(): Seq[Region] = {
    currentRegionsPriorityQueue.toSeq
  }

  /**
   * Advance to the specified locus, which must be greater than the current locus. After calling this, the
   * [[currentRegions]] method will give the overlapping regions at the new locus.
   *
   * @param locus Locus to advance to.
   * @return The *new regions* that were added as a result of this call. Note that this is not the full set of regions
   *         in the window: you must examine [[currentRegions]] for that.
   */
  def setCurrentLocus(locus: Long): Seq[Region] = {
    assume(locus >= currentLocus, "Pileup window can only move forward in locus")
    currentLocus = locus

    // Remove regions that are no longer in the window.
    // Note that the end of a region is exclusive, so e.g. if halfWindowSize=0 and head.end=locus, we do want to drop
    // that read.
    while (!currentRegionsPriorityQueue.isEmpty && currentRegionsPriorityQueue.head.end <= locus - halfWindowSize) {
      val dropped = currentRegionsPriorityQueue.dequeue()
      assert(!dropped.overlapsLocus(locus, halfWindowSize))
    }

    newRegions = if (sortedRegions.isEmpty) {
      Seq.empty
    } else {
      // Build up a list of new regions that are now in the window.
      // Note that the start of a region is inclusive, so e.g. if halfWindowSize=0 and head.start=locus we want to
      // include it.
      val newRegionsBuilder = mutable.ArrayBuffer.newBuilder[Region]
      while (sortedRegions.nonEmpty && sortedRegions.head.start <= locus + halfWindowSize) {
        val region = sortedRegions.next()
        if (region.overlapsLocus(locus, halfWindowSize)) newRegionsBuilder += region
      }
      newRegionsBuilder.result
    }
    currentRegionsPriorityQueue.enqueue(newRegions: _*)
    newRegions // We return the newly added regions.
  }

  /**
   * Find the next locus for which calling setCurrentLocus(locus) would result in a non-empty set of regions in
   * currentRegions.
   *
   * @return Some(locus) if such a locus exists, otherwise None
   */
  def nextLocusWithRegions(): Option[Long] = {
    if (currentRegionsPriorityQueue.exists(_.overlapsLocus(currentLocus + 1, halfWindowSize))) {
      Some(currentLocus + 1)
    } else if (sortedRegions.hasNext) {
      val result = math.max(0, sortedRegions.head.start - halfWindowSize)
      assert(result > currentLocus)
      Some(result)
    } else {
      None
    }
  }
}
object SlidingWindow {
  /**
   * Advance one or more sliding windows to the next locus in an iterator, optionally skipping loci for which no windows
   * have regions.
   *
   * If this function returns Some(locus), then all the windows are now positioned at locus. If it returns None, then
   * there are no more loci to process (we have reached the end of the loci iterator), and the positions of the sliding
   * windows are undefined.
   *
   * This function takes an iterator of loci instead of a single locus to enable it to efficiently skip entire ranges of
   * loci when skipEmpty=true.
   *
   * @param windows SlidingWindow instances to be advanced
   * @param loci iterator over loci. This function will advance the iterator at least once, and possibly many times
   *             if skipEmpty is true.
   * @param skipEmpty whether to skip over loci for which no window has any regions
   * @return Some(locus) if there was another locus left to process, otherwise None
   */
  def advanceMultipleWindows[Region <: HasReferenceRegion](windows: PerSample[SlidingWindow[Region]],
                                                           loci: LociSet.SingleContig.Iterator,
                                                           skipEmpty: Boolean = true): Option[Long] = {
    if (skipEmpty) {
      while (loci.hasNext) {
        val nextNonEmptyLocus = windows.flatMap(_.nextLocusWithRegions).reduceOption(_ min _)
        if (nextNonEmptyLocus.isEmpty) {
          // Our windows are out of regions. We're done.
          return None
        } else if (nextNonEmptyLocus.get <= loci.head) {
          // The next locus with regions is at or before the next locus in the iterator.
          // We advance to the next locus in the iterator, and check if the resulting windows are all empty.
          val nextLocus = loci.next()
          windows.foreach(_.setCurrentLocus(nextLocus))

          // Windows may still be empty here, because the next locus with regions may have been before the next locus,
          // and now we just fast-forwarded past the regions into an empty area of the genome.
          // If any window is nonempty, we're done. If not, we continue looping.
          if (windows.exists(_.currentRegions.nonEmpty)) {
            return Some(nextLocus)
          }
        } else {
          // The next locus with regions is beyond the next locus in the iterator.
          // We skip the iterator forward so in the next iteration of the loop, we'll be at the locus with regions.
          loci.skipTo(nextNonEmptyLocus.get)
        }
      }
      // No more loci remaining in the iterator. We're done.
      return None
    } else if (loci.hasNext) {
      // Not skipping empty, and we have another locus in the iterator to go to.
      val nextLocus = loci.next()
      windows.foreach(_.setCurrentLocus(nextLocus))
      return Some(nextLocus)
    } else {
      // We are out of loci.
      return None
    }
  }
}
