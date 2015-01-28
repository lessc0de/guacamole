package org.hammerlab.guacamole.alignment

import org.hammerlab.guacamole.Bases

case class Kmer(bases: Seq[Byte]) {
  def head = bases.head

  def tail = bases.tail

  def lastBase = bases.last

  def possibleNext: Seq[Kmer] = {
    Seq(
      Kmer(bases.tail :+ Bases.A),
      Kmer(bases.tail :+ Bases.T),
      Kmer(bases.tail :+ Bases.C),
      Kmer(bases.tail :+ Bases.G)
    )
  }

  override def toString: String = Bases.basesToString(bases)
}

object Kmer {
  def apply(seq: String): Kmer = {
    Kmer(Bases.stringToBases(seq).toArray)
  }
}
