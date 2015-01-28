package org.hammerlab.guacamole.alignment

import org.hammerlab.guacamole.Bases
import org.scalatest.{FunSuite, Matchers}

class PackedKmerSuite extends FunSuite with Matchers {

  test("short kmer packing: all A") {
    val allA = Bases.stringToBases("AAAA")
    PackedKmer(allA) should be(PackedKmer(0))
  }

  test("short kmer packing: all T") {

    val allT = Bases.stringToBases("TTTT")
    PackedKmer(allT) should be(PackedKmer(2 + 8 + 32 + 128))
  }

  test("short kmer packing: alternating A and T") {
    val altAT = Bases.stringToBases("ATAT")
    PackedKmer(altAT) should be(PackedKmer(2 + 0 + 32 + 0))
  }

  test("kmer last base") {
    val tcga = PackedKmer("TCGA")
    tcga.lastBaseAsNucleotide should be(Bases.A)

    val acgt = PackedKmer("ACGT")
    acgt.lastBaseAsNucleotide should be(Bases.T)
  }

  test("kmer tail") {
    val tcga = PackedKmer("TCGA")
    tcga.tail(3) should be(PackedKmer("CGA"))

    val agct = PackedKmer("AGCT")
    agct.tail(3) should be(PackedKmer("GCT"))
  }

  test("kmer head") {
    val tcga = PackedKmer("TCGA")
    tcga.head should be(PackedKmer("TCG"))

    val acgt = PackedKmer("ACGT")
    acgt.head should be(PackedKmer("ACG"))
  }

  test("kmer add") {
    val allA = PackedKmer("AAAA")
    val changeTailToT = allA + 2
    changeTailToT should be (PackedKmer("AAAT"))
  }

  test("kmer multiply") {
    val aat = PackedKmer("AAT")
    val shifted = aat * 4
    shifted should be (PackedKmer("AATA"))
  }

  test("kmer advance") {
    val aat = PackedKmer("AAT")
    val shifted = aat.advance(3, Bases.T)
    shifted should be (PackedKmer("ATT"))
    shifted.advance(3, Bases.C) should be(PackedKmer("TTC"))
  }

  test("multiple kmers from long string") {
    val seq = Bases.stringToBases("ATATCCGG")
    val allKmers = PackedKmer(seq, 3)

    allKmers.length should be (6)
    allKmers(0) should be (PackedKmer("ATA"))
    allKmers(1) should be (PackedKmer("TAT"))
    allKmers(2) should be (PackedKmer("ATC"))
    allKmers(3) should be (PackedKmer("TCC"))
    allKmers(4) should be (PackedKmer("CCG"))
    allKmers(5) should be (PackedKmer("CGG"))

  }

  test("max length of Kmer") {
    val seq = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
    println(seq.length)
    val long = PackedKmer(seq)
    long.packedKmer should not be (0)
    assert(long.packedKmer > 0)
  }

}
