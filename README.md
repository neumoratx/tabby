# Tabby

Tabby is a project heavily based on htslib's `Tabix/bgzip` tools for the purpose of more general but still efficient querying of tabular data on the command line.

The idea is to extend the already excellent core capabilities of the Tabix indexing scheme to allow for more flexible indexing and filtering (querying) of multiple columns all in one call to Tabby.   While Tabix only did chromosome + position indexing, Tabby switches that to general indexing + filtering operations that can be used across non-genomic datatypes.  

A primary design goal is putting the least amount of burden on the user as possible for small-to-medium sized, row-major data that can fit on one or a few machines (i.e. doesn't truly need a cluster to compute over).

Thus, while other tools may boast faster retrieval, better compression,  and/or more extensive filtering/querying expressions (e.g. SQL-support), they do so at the cost of conversion to a non-original, non-flatfile format without default compression (e.g. SQLite/PostgreSQL), non-row-major focus with a more extensive set of dependencies/code (e.g. Parquet), and/or external service subscriptions/opaque implementations (Snowflake/BigQuery/RedShift).

NOTE: an LLM did most/all code changes, though under careful human scrutiny and constantly human-guided over many prompts.

## Why?
Some might ask what value is there is in doing this.  

Let us start with the following observations about the beneficial features of Tabix as it's currently implemented:
* Fast random-access, row-based retrieval
* Compression by default
* Support for retrieval on remote files on HTTP/HTTPS/S3/GCS stores
* File format is nearly original (tab-delimited text files compressed using a gzip compatible method)

The key observation here is that we get to essentially keep our original tabular file format, but now have we can do random access queries against a compressed version of it that's stored on an object store (or other remote site) w/o having to download and decompress the full thing ahead of time and without having to write any code to do this.

However, this use case is limited to genomic region-based querying, period.

Might there be other uses for the above benefits, while not being tied to a genomics-region only context? Indeed there are.

The following 3 main features above/beyond Tabix are detailed below. 

Please keep in mind that you can use any or all of these, they don't require each other to work.  

Nor do you need to change an existing Tabix index (`.tbi` file) to use the regular region querying with additional-column filtering.

## Single Column Indexing
Tabby by default supports indexing on a single column instead of the chromosome + start + end normally required Tabix (`-p seqonly` instead of `-s 1 -b 2 -e 3`).  

This was a simple modification to leverage the existing b-tree index used for the chromosome and apply it instead to any column that you've already sorted on, no longer needing a position column as well.

While simple to make, this change allows a lot more flexibility w/o resorting to crude hacks where a "stub" column for a position is added to the file and the chromosome column is overriden to your desired column to index on.

One common use case (at least of the author) is doing this for gene names or gene IDs, instead of their coordinates.

## Multi-Column Filtering
Tabby also adds support for one or more additional filters on the columns in the file, supporting both numeric (equal `==`, not-equal `!=`, greater-than-equal `>=`, lesser-than-equal `<=`) and string operations (equal `==`, not-equal `!=`, and basic regex: `~=`  and/or `!~`).  

Both AND (`-F`) and OR (`-O`) are supported when specifying a specific filter.

Columns can be specified by either their *exact* name (case sensitive) or by column index (starting at 0 from the left).

Column names if used in the filters *must not* contain any of the 2-character operators: `==`, `!=`, `>=`, `<=`, `~=`, nor `!~`.

## Secondary Index Across All Numeric Columns
An additional important, but more "under-the-hood" feature that Tabby has beyond Tabix, is that while the above filtering will work on a normal bgzipped file with normal tabix index, Tabby now supports the generation of a secondary index (`sidx`) which stores additional information across all the numeric columns (not applicable to string columns).  This secondary index specifically tracks the min/max values of each numeric column across each compressed block in the block-gzip.  

Behond providing potentially useful information about the structure of the original file + index, when a numeric filter is applied and the secondary index exists, Tabby can leverage this to drop blocks (and potentially chunks of the file in a remote file context, *if* the chunk has *no* blocks which match any of the query/filters) from processing (decompressing + reading line-by-line).

This can potentially speed up processing by further pruning out whole blocks (or even whole chunks) from consideration*.

All that to say, Tabby is *not* making any major performance enhancement(s) to the core htslib code. 

That said, we think having multi-column filtering support inline with the single column index can save performance when otherwise you might need to combine a `tabix` call with one or more piped `grep`s to get the same behavior.

*in a remote file context the performance gains are more minimal, if any, due to how `htslib` handles retrieval of blocks.  Retrieval in this context is at thepotentially much larger chunk level to save costly remote calls, rather than at the more fine-grained block level.  Since chunks could have at least one block that matches one or more of the filters, it's more likely you'll download the full chunk.  Decompression can still be skipped for individual blocks though.

### Building HTSlib (and by extension Tabby)

See [INSTALL](INSTALL) for complete details.
[Release tarballs][download] contain generated files that have not been
committed to this repository, so building the code from a Git repository
requires extra steps:

```sh
autoreconf -i  # Build the configure script and install files it uses
./configure    # Optional but recommended, for choosing extra functionality
make
make install
```

[download]: http://www.htslib.org/download/

### Citing

Please cite this paper when using HTSlib for your publications.

> HTSlib: C library for reading/writing high-throughput sequencing data </br>
> James K Bonfield, John Marshall, Petr Danecek, Heng Li, Valeriu Ohan, Andrew Whitwham, Thomas Keane, Robert M Davies </br>
> _GigaScience_, Volume 10, Issue 2, February 2021, giab007, https://doi.org/10.1093/gigascience/giab007

```
@article{10.1093/gigascience/giab007,
    author = {Bonfield, James K and Marshall, John and Danecek, Petr and Li, Heng and Ohan, Valeriu and Whitwham, Andrew and Keane, Thomas and Davies, Robert M},
    title = "{HTSlib: C library for reading/writing high-throughput sequencing data}",
    journal = {GigaScience},
    volume = {10},
    number = {2},
    year = {2021},
    month = {02},
    abstract = "{Since the original publication of the VCF and SAM formats, an explosion of software tools have been created to process these data files. To facilitate this a library was produced out of the original SAMtools implementation, with a focus on performance and robustness. The file formats themselves have become international standards under the jurisdiction of the Global Alliance for Genomics and Health.We present a software library for providing programmatic access to sequencing alignment and variant formats. It was born out of the widely used SAMtools and BCFtools applications. Considerable improvements have been made to the original code plus many new features including newer access protocols, the addition of the CRAM file format, better indexing and iterators, and better use of threading.Since the original Samtools release, performance has been considerably improved, with a BAM read-write loop running 5 times faster and BAM to SAM conversion 13 times faster (both using 16 threads, compared to Samtools 0.1.19). Widespread adoption has seen HTSlib downloaded \\&gt;1 million times from GitHub and conda. The C library has been used directly by an estimated 900 GitHub projects and has been incorporated into Perl, Python, Rust, and R, significantly expanding the number of uses via other languages. HTSlib is open source and is freely available from htslib.org under MIT/BSD license.}",
    issn = {2047-217X},
    doi = {10.1093/gigascience/giab007},
    url = {https://doi.org/10.1093/gigascience/giab007},
    note = {giab007},
    eprint = {https://academic.oup.com/gigascience/article-pdf/10/2/giab007/36332285/giab007.pdf},
}
```

### Support

If you have found a bug or would like a new feature, please report the same in the GitHub [HTSlib](https://github.com/samtools/htslib/issues) issue tracker.

For any security related issue, please send a mail to [samtools@sanger.ac.uk](mailto:samtools@sanger.ac.uk) instead of reporting in the GitHub issue tracker. 
