package org.hammerlab.guacamole.alignment

import org.hammerlab.guacamole.Bases

import scala.collection.mutable.ArrayBuffer

case class Kmer(packedKmer: Long) extends AnyVal {

  def lastBaseAsNucleotide: Byte = {
    Kmer.BASE_VALUES_REVERSE(lastBase)
  }

  def lastBase: Int = {
    (packedKmer % Kmer.NUM_BASES).toInt
  }

  // All but the last base: (n - 1) length prefix
  def head: Kmer = {
    Kmer((packedKmer - lastBase) / Kmer.NUM_BASES)
  }

  // All but the first base: (size - 1) length suffix
  def tail(tailSize: Int): Kmer = {
    Kmer(packedKmer % Kmer.ROOT_BASES(tailSize))
  }

  def *(x: Int): Kmer = {
    Kmer(packedKmer * x)
  }

  def +(x: Int): Kmer = {
    Kmer(packedKmer + x)
  }

  def possibleNext(size: Int): Iterable[Kmer] = {
    val tailKmer = tail(size - 1) * Kmer.NUM_BASES
    Kmer.BASE_VALUES.values.map(b => tailKmer + b)
  }

  def advance(size: Int, nextBase: Byte): Kmer = {
    (tail(size - 1) * Kmer.NUM_BASES) + Kmer.BASE_VALUES(nextBase)
  }

  def append(size: Int, nextBase: Byte): Kmer = {
    Kmer((packedKmer * Kmer.NUM_BASES) + Kmer.BASE_VALUES(nextBase))
  }

  def toBases(size: Int): Seq[Byte] = {
    val result = ArrayBuffer.newBuilder[Byte]
    var packedVal = packedKmer

    var i = 0
    while (i < size) {
      result += Kmer.BASE_VALUES_REVERSE((packedVal % Kmer.NUM_BASES).toInt)
      packedVal /= Kmer.NUM_BASES
      i += 1
    }
    result.result().reverse
  }

  def toString(size: Int): String = {
    Bases.basesToString(toBases(size))
  }
}

object Kmer {

  val BASE_VALUES = Map(
    (Bases.A, 0),
    (Bases.C, 1),
    (Bases.T, 2),
    (Bases.G, 3)
  )

  val BASE_VALUES_REVERSE = Array(Bases.A, Bases.C, Bases.T, Bases.G)
  val NUM_BASES = BASE_VALUES.size
  val ROOT_BASES = (0 until 31).map(math.pow(NUM_BASES, _).toInt).toArray

  def apply(seq: Seq[Byte]): Kmer = {
    Kmer(getKmerPackedVal(seq))
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

  def apply(seq: String): Kmer = {
    Kmer(Bases.stringToBases(seq))
  }

  def apply(seq: Seq[Byte], kmerSize: Int): Seq[Kmer] = {

    val numKmers = seq.size - kmerSize + 1
    var lastKmer = Kmer(seq.take(kmerSize))
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

  def kmerPathToSeq(kmerPath: Seq[Kmer], kmerSize: Int): Seq[Byte] = {
    val result = ArrayBuffer.newBuilder[Byte]
    result ++= kmerPath.head.toBases(kmerSize)

  }

  def kmerPathToString(kmerPath: Seq[Kmer]): String] = {
    Bases.basesToString(kmerPathToSeq(kmerPath))
  }
}
