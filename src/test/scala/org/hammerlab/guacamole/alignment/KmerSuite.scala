package org.hammerlab.guacamole.alignment

import org.hammerlab.guacamole.Bases
import org.scalatest.{FunSuite, Matchers}

class KmerSuite extends FunSuite with Matchers {

  test("short kmer packing: all A") {
    val allA = Bases.stringToBases("AAAA")
    Kmer(allA) should be(Kmer(0))
  }

  test("short kmer packing: all T") {

    val allT = Bases.stringToBases("TTTT")
    Kmer(allT) should be(Kmer(2 + 8 + 32 + 128))
  }

  test("short kmer packing: alternating A and T") {
    val altAT = Bases.stringToBases("ATAT")
    Kmer(altAT) should be(Kmer(2 + 0 + 32 + 0))
  }

  test("kmer last base") {
    val tcga = Kmer("TCGA")
    tcga.lastBaseAsNucleotide should be(Bases.A)

    val acgt = Kmer("ACGT")
    acgt.lastBaseAsNucleotide should be(Bases.T)
  }

  test("kmer tail") {
    val tcga = Kmer("TCGA")
    tcga.tail(3) should be(Kmer("CGA"))

    val agct = Kmer("AGCT")
    agct.tail(3) should be(Kmer("GCT"))
  }

  test("kmer head") {
    val tcga = Kmer("TCGA")
    tcga.head should be(Kmer("TCG"))

    val acgt = Kmer("ACGT")
    acgt.head should be(Kmer("ACG"))
  }

  test("kmer add") {
    val allA = Kmer("AAAA")
    val changeTailToT = allA + 2
    changeTailToT should be (Kmer("AAAT"))
  }

  test("kmer multiply") {
    val aat = Kmer("AAT")
    val shifted = aat * 4
    shifted should be (Kmer("AATA"))
  }

  test("kmer advance") {
    val aat = Kmer("AAT")
    val shifted = aat.advance(3, Bases.T)
    shifted should be (Kmer("ATT"))
    shifted.advance(3, Bases.C) should be(Kmer("TTC"))
  }

  test("multiple kmers from long string") {
    val seq = Bases.stringToBases("ATATCCGG")
    val allKmers = Kmer(seq, 3)

    allKmers.length should be (6)
    allKmers(0) should be (Kmer("ATA"))
    allKmers(1) should be (Kmer("TAT"))
    allKmers(2) should be (Kmer("ATC"))
    allKmers(3) should be (Kmer("TCC"))
    allKmers(4) should be (Kmer("CCG"))
    allKmers(5) should be (Kmer("CGG"))

  }

  test("max length of Kmer") {
    val seq = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
    println(seq.length)
    val long = Kmer(seq)
    long.packedKmer should not be (0)
    assert(long.packedKmer > 0)
  }

}
