#!/usr/bin/env python
"""
discordex.py

Merges all UNCID_*.sorted_genome_alignments.bam.discord.tsv.gz files containing
instances of discordant alignments rounded to nearest 100th coordinate.
Considers only primary alignments.

Requires directory containing output of omfgene_wrapper from CGC. See README.md
for more information on how to reproduce this output.

Also requires TCGA.tsv, which is obtained from
http://duffel.rail.bio/recount/TCGA/TCGA.tsv .

Writes two files:
samples.tsv, with tab-separated columns
1. sample index
2. GDC case ID
3. GDC file ID for tar
4. BAM filename

discordex.tsv, with tab separated columns
1. mate 1 chrom
2. mate 1 alignment position (1-based, rounded to nearest 100)
3. mate 2 chrom
4. mate 2 alignment position (1-based, rounded to nearest 100)
5. comma-separated list of sample indexes
6. comma-separated list of corresponding coverages

We ran

pypy discordex.py --tcga-metadata /path/to/TCGA.tsv
                  --omfgene-output /path/to/omfgene/output
                  --output /path/to/output/dir
"""
import glob
import os
import sys
import errno
import tempfile
import gzip
import contextlib
import shutil
import atexit
import subprocess

@contextlib.contextmanager
def xopen(gzipped, *args):
    """ Passes args on to the appropriate opener, gzip or regular.

        In compressed mode, functionality almost mimics gzip.open,
        but uses gzip at command line.

        As of PyPy 2.5, gzip.py appears to leak memory when writing to
        a file object created with gzip.open().

        gzipped: True iff gzip.open() should be used to open rather than
            open(); False iff open() should be used; None if input should be
            read and guessed; '-' if writing to stdout
        *args: unnamed arguments to pass

        Yield value: file object
    """
    import sys
    if gzipped == '-':
        fh = sys.stdout
    else:
        if not args:
            raise IOError('Must provide filename')
        if gzipped is None:
            with open(args[0], 'rb') as binary_input_stream:
                # Check for magic number
                if binary_input_stream.read(2) == '\x1f\x8b':
                    gzipped = True
                else:
                    gzipped = False
        if gzipped:
            try:
                mode = args[1]
            except IndexError:
                mode = 'rb'
            if 'r' in mode:
                # Be forgiving of gzips that end unexpectedly
                old_read_eof = gzip.GzipFile._read_eof
                gzip.GzipFile._read_eof = lambda *args, **kwargs: None
                fh = gzip.open(*args)
            elif 'w' in mode or 'a' in mode:
                try:
                    compresslevel = int(args[2])
                except IndexError:
                    compresslevel = 9
                if 'w' in mode:
                    output_stream = open(args[0], 'wb')
                else:
                    output_stream = open(args[0], 'ab')
                gzip_process = subprocess.Popen(['gzip',
                                                    '-%d' % compresslevel],
                                                    bufsize=-1,
                                                    stdin=subprocess.PIPE,
                                                    stdout=output_stream)
                fh = gzip_process.stdin
            else:
                raise IOError('Mode ' + mode + ' not supported')
        else:
            fh = open(*args)
    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()
        if 'gzip_process' in locals():
            gzip_process.wait()
        if 'output_stream' in locals():
            output_stream.close()
        if 'old_read_eof' in locals():
            gzip.GzipFile._read_eof = old_read_eof

class xstream(object):
    """ Permits Pythonic iteration through partitioned/sorted input streams.

        All iterators are implemented as generators. Could have subclassed
        itertools.groupby here; however, implementation of itertools.groupby
        may change from version to version of Python. Implementation is thus
        just based on itertools.groupby from
        https://docs.python.org/2/library/itertools.html .

        Usage: for key, xpartition in xstream(hadoop_stream):
                   for value in xpartition:
                        <code goes here>

        Each of key and value above is a tuple of strings.

        Properties
        -------------
        key: key tuple that denotes current partition; this is an attribute
            of both an xstream. None when no lines have been read yet.
        value: tuple that denotes current value. None when no lines have been
            read yet.

        Init vars
        -------------
        input_stream: where to find input lines
        key_fields: the first "key_fields" fields from an input line are
            considered the key denoting a partition
        separator: delimiter separating fields from each input line
        skip_duplicates: skip any duplicate lines that may follow a line
    """
    @staticmethod
    def stream_iterator(
            input_stream,
            separator='\t',
            skip_duplicates=False
        ):
        if skip_duplicates:
            for line, _ in groupby(input_stream):
                yield tuple(line.strip().split(separator))
        else:
            for line in input_stream:
                yield tuple(line.strip().split(separator))

    def __init__(
            self, 
            input_stream,
            key_fields=1,
            separator='\t',
            skip_duplicates=False
        ):
        self._key_fields = key_fields
        self.it = self.stream_iterator(
                        input_stream,
                        separator=separator,
                        skip_duplicates=skip_duplicates
                    )
        self.tgtkey = self.currkey = self.currvalue = object()

    def __iter__(self):
        return self

    def next(self):
        while self.currkey == self.tgtkey:
            self.currvalue = next(self.it)    # Exit on StopIteration
            self.currkey = self.currvalue[:self._key_fields]
        self.tgtkey = self.currkey
        return self.currkey, self._grouper(self.tgtkey)

    def _grouper(self, tgtkey):
        while self.currkey == tgtkey:
            yield self.currvalue[self._key_fields:]
            self.currvalue = next(self.it)    # Exit on StopIteration
            self.currkey = self.currvalue[:self._key_fields]

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--tcga-metadata', type=str, required=True,
        help='path to TCGA metadata file, which can be downloaded at '
             'http://duffel.rail.bio/recount/TCGA/TCGA.tsv')
    parser.add_argument('--omfgene-output', type=str, required=True,
        help='directory containing omfgene tsv.gzs')
    parser.add_argument('--output', type=str, required=True,
        help='output directory')
    parser.add_argument('--temp-dir', type=str, required=False,
        default=None,
        help='where to store temporary files')
    args = parser.parse_args()
    '''Consider only UNCID*.bam.discord.tsv.gz's; others correspond to 
    different alignment protocols.'''
    omfgene_files = glob.glob(
            os.path.join(args.omfgene_output,
                         'UNCID*.sorted_genome_alignments.bam.discord.tsv.gz')
        )
    bam_id_to_gdc_uuid = {}
    with open(args.tcga_metadata) as metadata_stream:
        metadata_stream.readline()
        for line in metadata_stream:
            tokens = line.strip().split('\t')
            try:
                bam_id_to_gdc_uuid[tokens[24].split('.')[1]] = (
                        tokens[21], tokens[22]
                    )
            except IndexError:
                print >>sys.stderr, (
                        'Warning: unexpected number of columns in line "{}"'
                    ).format(line.strip())
    try:
        os.makedirs(args.output)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    omfgene_file_to_index = {}
    with open(os.path.join(args.output, 'samples.tsv'), 'w') as index_stream:
        print >>index_stream, '\t'.join(
                ['sample index', 'case uuid', 'tar.gz uuid', 'bam filename']
            )
        i = 0
        for omfgene_file in omfgene_files:
            split_name = omfgene_file.split('.')
            try:
                id_info = bam_id_to_gdc_uuid[split_name[1]]
            except KeyError:
                id_info = ('NA', 'NA')
            omfgene_file_to_index[omfgene_file] = i
            bam = '.'.join(split_name[:-3])
            print >>index_stream, '\t'.join(
                    map(str,
                            [i, id_info[0], id_info[1], os.path.basename(bam)]
                        )
                )
            i += 1
    # Merge
    temp_dir = tempfile.mkdtemp(dir=args.temp_dir)
    atexit.register(shutil.rmtree, temp_dir)
    with xopen(True, os.path.join(
                temp_dir, 'to_merge.tsv.gz',
            ), 'w') as merge_stream:
        for i, omfgene_file in enumerate(omfgene_file_to_index):
            print 'Writing file {} to merged index.\r'.format(i),
            sys.stdout.flush()
            index = str(omfgene_file_to_index[omfgene_file])
            with gzip.open(omfgene_file) as omfgene_stream:
                for line in omfgene_stream:
                    line = line.strip()
                    if not line.endswith('0'): continue
                    tokens = line.split(' ')
                    tokens = tokens[1].split('\t')[:-1] + [tokens[0], index]
                    print >>merge_stream, '\t'.join(tokens)

    sort_process = subprocess.Popen(
        ('set -exo pipefail; '
         'cd {}; gzip -cd {} | '
         'sort -k1,1 -k2,2n -k3,3 -k4,4n -k6,6n | gzip >{}').format(
                temp_dir, 'to_merge.tsv.gz',
                'merged.tsv.gz'
            ), executable='/bin/bash', shell=True)
    sort_process.wait()
    with xopen(
            None, os.path.join(temp_dir, 'merged.tsv.gz')
        ) as merged_stream, xopen(True,
                os.path.join(args.output, 'discordex.v1.hg19.tsv.gz'
            ), 'w') as discordex_stream:
        for key, xpartition in xstream(merged_stream, 4):
            samples_and_coverages = zip(*list(xpartition))[::-1]
            print >>discordex_stream, '\t'.join(
                    key
                    + (','.join(samples_and_coverages[0]),
                       ','.join(samples_and_coverages[1]))
                )
