guacamole
=========
[![Build Status](https://travis-ci.org/hammerlab/guacamole.svg?branch=master)](https://travis-ci.org/hammerlab/guacamole)

Guacamole is a framework for variant calling, i.e. identifying DNA mutations
from [Next Generation Sequencing](http://en.wikipedia.org/wiki/DNA_sequencing)
data. It currently includes a toy germline (non-cancer) variant caller as well
as a somatic variant caller for finding cancer mutations.  Most development
effort has gone into the somatic caller so far.

The emphasis is on a readable codebase that can be readily understood and
adapted for experimentation.

Guacamole is written in Scala using the [Apache
Spark](http://spark.apache.org/) engine for distributed processing. It can run
on a single computer or on a Hadoop cluster.

Guacamole supports reading aligned reads as:
 * [BAM or SAM files](http://genomicsandhealth.org/our-work/work-products/file-formats-sambam)
 * [Parquet files](http://parquet.incubator.apache.org/) of [bdg-formats](https://github.com/bigdatagenomics/bdg-formats) [AlignmentRecords](https://github.com/bigdatagenomics/bdg-formats/blob/master/src/main/resources/avro/bdg.avdl#L60)
 
It can write the called genotypes as:
 * [VCF files](http://genomicsandhealth.org/our-work/work-products/file-formats-vcfbcf)
 * Parquet files of bdg-formats [Genotypes](https://github.com/bigdatagenomics/bdg-formats/blob/master/src/main/resources/avro/bdg.avdl#L547).

Guacamole uses ideas and some functionality from
[ADAM](https://github.com/bigdatagenomics/adam). It also takes inspiration from
the [Avocado](https://github.com/bigdatagenomics/avocado) project.

For hacking Guacamole, see our [code docs](http://www.hammerlab.org/guacamole/docs/#org.hammerlab.guacamole.package).

# Running Guacamole on a Single Node

Guacamole requires [Apache Maven](http://maven.apache.org/).

Build:

```
mvn package
```

This will build a guacamole JAR file in the `target` directory. A script is
included to run it:

```
scripts/guacamole germline-threshold \
	-reads src/test/resources/chrM.sorted.bam \
	-out /tmp/result.vcf
```

This creates a *directory* called `/tmp/result.vcf`. The actual VCF file is in
`/tmp/result.vcf/part-r-00000`. You'll always get one part file in the output
directory.

Try 
```
scripts/guacamole -h
```
for a list of implemented variant callers, or

```
scripts/guacamole <caller> -h
```
for help on a particular variant caller.

# Running Guacamole on a Hadoop Cluster

See Guacamole's
[pom.xml](https://github.com/hammerlab/guacamole/blob/master/pom.xml) file for
the versions of Hadoop and Spark that Guacamole expects to find on your
cluster.

Here is an example command to get started using Guacamole in Spark's yarn
cluster mode. You'll probably have to modify it for your environment. 

```
spark-submit \
	--master yarn \
	--deploy-mode cluster \
	--driver-java-options -Dlog4j.configuration=/path/to/guacamole/scripts/log4j.properties \
	--executor-memory 4g \
	--driver-memory 10g \
	--num-executors 1000 \
	--executor-cores 1 \
	--class org.hammerlab.guacamole.Guacamole \
	--verbose \
	/path/to/target/guacamole-0.0.1.jar \
	germline-threshold \
        -reads hdfs:///path/to/reads.bam \
        -out hdfs:///path/to/result.vcf \
	-spark_master yarn-cluster
```

# Is this ready for production use?

Not currently. Everything here is experimental. Please use a standard tool if
you need accurate variant calls.


# License

Guacamole is released under an [Apache 2.0 license](LICENSE.txt).

YourKit is kindly supporting this open source project with its full-featured Java Profiler.
YourKit, LLC is the creator of innovative and intelligent tools for profiling
[Java](http://www.yourkit.com/java/profiler/index.jsp) and [.NET](http://www.yourkit.com/.net/profiler/index.jsp) applications.
![image](https://cloud.githubusercontent.com/assets/455755/4988560/97757f12-6935-11e4-9270-f5fc42f9b585.png)
