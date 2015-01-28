package org.hammerlab.guacamole.alignment

import org.hammerlab.guacamole.Bases

import scala.collection.mutable.ArrayBuffer

/**
 * Kmer packed into a 64-bit long with 2bits/base (supports up to 31 length Kmer)
 * @param packedKmer
 */
case class PackedKmer(packedKmer: Long) extends AnyVal {

  /**
   * Last base in the kmer (as Bases.A/T/C/G)
   */
  def lastBaseAsNucleotide: Byte = {
    PackedKmer.BASE_VALUES_REVERSE(lastBase)
  }

  /**
   *  Last base in the kmer (coded as 0/1/2/3)
   * @return
   */
  def lastBase: Int = {
    (packedKmer % PackedKmer.NUM_BASES).toInt
  }

  // All but the last base: (n - 1) length prefix
  def head: PackedKmer = {
    PackedKmer((packedKmer - lastBase) / PackedKmer.NUM_BASES)
  }

  // All but the first base: (size - 1) length suffix
  def tail(tailSize: Int): PackedKmer = {
    PackedKmer(packedKmer % PackedKmer.ROOT_BASES(tailSize))
  }

  def *(x: Int): PackedKmer = {
    PackedKmer(packedKmer * x)
  }

  def +(x: Int): PackedKmer = {
    PackedKmer(packedKmer + x)
  }

  def possibleNext(size: Int): Iterable[PackedKmer] = {
    val tailKmer = tail(size - 1) * PackedKmer.NUM_BASES
    PackedKmer.BASE_VALUES.values.map(b => tailKmer + b)
  }

  def advance(size: Int, nextBase: Byte): PackedKmer = {
    (tail(size - 1) * PackedKmer.NUM_BASES) + PackedKmer.BASE_VALUES(nextBase)
  }

  def append(size: Int, nextBase: Byte): PackedKmer = {
    PackedKmer((packedKmer * PackedKmer.NUM_BASES) + PackedKmer.BASE_VALUES(nextBase))
  }

  def toBases(size: Int): Seq[Byte] = {
    val result = ArrayBuffer.newBuilder[Byte]
    var packedVal = packedKmer

    var i = 0
    while (i < size) {
      result += PackedKmer.BASE_VALUES_REVERSE((packedVal % PackedKmer.NUM_BASES).toInt)
      packedVal /= PackedKmer.NUM_BASES
      i += 1
    }
    result.result().reverse
  }

  def toString(size: Int): String = {
    Bases.basesToString(toBases(size))
  }
}

object PackedKmer {

  val BASE_VALUES = Map(
    (Bases.A, 0),
    (Bases.C, 1),
    (Bases.T, 2),
    (Bases.G, 3)
  )

  val BASE_VALUES_REVERSE = Array(Bases.A, Bases.C, Bases.T, Bases.G)
  val NUM_BASES = BASE_VALUES.size
  val ROOT_BASES = (0 until 31).map(math.pow(NUM_BASES, _).toInt).toArray

  def apply(seq: Seq[Byte]): PackedKmer = {
    PackedKmer(getKmerPackedVal(seq))
  }

  def getKmerPackedVal(seq: Seq[Byte]): Long = {
    assume(seq.size < 32, "Kmer exceeds max kmer length (31)")
    Bases.assertAllStandardBases(seq)

    var fourBase =  1L
    seq
      .reverse
      .zipWithIndex
      .map( {case (base, i) =>
      val next = BASE_VALUES(base) * fourBase
      fourBase *= NUM_BASES
      next
    }).sum
  }

  def apply(seq: String): PackedKmer = {
    PackedKmer(Bases.stringToBases(seq))
  }

  def apply(seq: Seq[Byte], kmerSize: Int): Seq[PackedKmer] = {

    val numKmers = seq.size - kmerSize + 1
    var lastKmer = PackedKmer(seq.take(kmerSize))
    val packedKmers = ArrayBuffer(lastKmer)
    var i = 1
    while (i < numKmers) {
      val nextKmer = lastKmer.advance(kmerSize, seq(kmerSize + i - 1))
      packedKmers += nextKmer
      lastKmer = nextKmer
      i += 1
    }
    packedKmers.toSeq
  }

  def kmerPathToSeq(kmerPath: Seq[PackedKmer], kmerSize: Int): Seq[Byte] = {
    val result = ArrayBuffer.newBuilder[Byte]
    result ++= kmerPath.head.toBases(kmerSize)
    result.result()
  }

  def kmerPathToString(kmerPath: Seq[PackedKmer], kmerSize: Int): String = {
    Bases.basesToString(kmerPathToSeq(kmerPath, kmerSize))
  }
}
