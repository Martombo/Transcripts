import features as ft
import subprocess as sp
import pysam as ps
import os


class Gtf:
    """utility functions to parse gtf"""

    def __init__(self, gtf_file, gene_biotypes=None, test=False):
        """
        :param gtf_file: path to gtf file
        :param gene_biotypes: list of gene biotypes to be considered
        :param test: for testing
        """
        if not test:
            assert isinstance(gtf_file, str)
            assert os.path.isfile(gtf_file)
        self.gtf_path = gtf_file
        self.gene_biotype = gene_biotypes
        self.genome = ft.Genome()

    def get_genome(self, fasta_path='', predict_noncoding_cds=False):
        """returns a Genome object, from an Ensembl gtf file"""
        self.genome.fasta_path = fasta_path
        with open(self.gtf_path) as f:
            for linea in f.readlines():
                line_dic = self._parse_line(linea)
                gene = self._get_gene(line_dic)
                if not line_dic or 'transcript_id' not in line_dic or \
                        (self.gene_biotype and line_dic['gene_biotype'] not in self.gene_biotype):
                    continue
                trans = self._get_transcript(line_dic, gene)
                if line_dic['type'] == 'start_codon':
                    trans.add_cds(ft.fix_order(line_dic['start'], line_dic['stop'], trans.strand)[0])
                elif line_dic['type'] == 'stop_codon':
                    trans.cds_stop = ft.fix_order(line_dic['start'], line_dic['stop'], trans.strand)[0]
                elif line_dic['type'] == 'exon':
                    self._set_exon(line_dic['start'], line_dic['stop'], trans, gene, int(line_dic['exon_number']), line_dic['exon_id'])
                elif line_dic['type'] == 'CDS':
                    trans.add_cds(line_dic['start'], line_dic['stop'])
        if predict_noncoding_cds:
            self.set_noncoding_cds()
        return self.genome

    def set_noncoding_cds(self):
        for trans in self.genome.transcripts:
            if not trans.cds_start:
                continue
            cds_start_exon = trans.exon_with(trans.cds_start)
            cds_start_exon.cds_start = trans.cds_start
            self._add_cds_from_exon(cds_start_exon)

    @staticmethod
    def _add_cds_from_exon(exon):
        for trans in exon.transcripts:
            if trans.cds_stop:
                continue
            trans.set_cds_from_start(exon.cds_start)

    @staticmethod
    def _set_exon(start, stop, trans, gene, num, id):
        if (start, stop) in gene.exons_dict:
            exon = gene.exons_dict[(start, stop)]
            exon.add_transcript(trans, num)
        else:
            exon = ft.Exon(id, num, trans, start, stop)
            gene.add_exon(exon, start, stop)

    @staticmethod
    def _get_transcript(attr, gene):
        """
        returns a Transcript from the Gene, if it's already present.
        otherwise a new Transcript
        """
        trans_id = attr['transcript_id']
        trans_biotype = attr['transcript_biotype']
        if trans_id in gene.trans_dict:
            return gene.trans_dict[trans_id]
        return ft.Transcript(trans_id, gene, biotype=trans_biotype)

    def _get_gene(self, dic):
        """
        returns a Gene from the Genome, if it's already present.
        otherwise a new Gene
        """
        if not dic:
            return
        chrom = self._get_chrom(dic['chr'])
        gene_id = dic['gene_id']
        if gene_id in chrom.genes_dict:
            return chrom.genes_dict[gene_id]
        gene = ft.Gene(gene_id, chrom, dic['gene_name'], dic['strand'])
        return gene

    def _get_chrom(self, chrom_name):
        """
        returns a Chromosome from the Genome, if it's already present.
        otherwise a new Chromosome
        """
        if chrom_name in self.genome.chroms_dict:
            return self.genome.chroms_dict[chrom_name]
        return ft.Chromosome(chrom_name, self.genome)

    def _parse_line(self, string):
        """parses a gtf line, with attributes (see self._parse_attributes())
        :rtype: dict
        :returns {chr, annot, type, start, ..., attr}"""
        splat = string.rstrip('\n').split('\t')
        if len(splat) < 8:
            return
        dic = dict(chr=splat[0], annot=splat[1], type=splat[2], start=int(splat[3]), stop=int(splat[4]),
                        score=splat[5], strand=splat[6], frame=splat[7])
        return self._add_attributes(dic, splat[8])

    @staticmethod
    def _add_attributes(dic, attrs):
        """parses the attributes in a dict
        :returns {'gene_id': 'ENSG000000123', 'gene_version' = '1', ...}
        """
        attrs_splat = attrs.split(';')
        for attr in attrs_splat:
            if not attr:
                continue
            attr = attr.lstrip(' ')
            attr_key = attr.split(' ')[0]
            attr_item = attr.split('"')[1]
            dic[attr_key] = attr_item
        return dic


class Bam:
    """utility functions to parse bam"""

    def __init__(self, path=None, reads_orientation='forward', test=False):
        """
        positions are 0-based
        :param reads_orientation: either 'forward', 'reverse' or 'mixed'
        """
        if not test:
            assert path
            assert path[-4:] == '.bam'
            assert os.path.isfile(path)
            if not os.path.isfile(path + '.bai'):
                p_index = sp.Popen(['samtools', 'index', path])
                p_index.communicate()
            assert os.path.isfile(path + '.bai')
            assert reads_orientation in ['forward', 'reverse', 'mixed']
            self.pysam = ps.AlignmentFile(path, 'rb')
        self.reads_orientation = reads_orientation

    def get_read_starts(self, chrom, start, stop, strand='', min_qual=40, read_len=[]):
        """
        collects the read start sites from the specified interval
        :param chrom:  str chromosome name
        :param start: int start 0-based
        :param stop: int stop 0-based
        :param strand: '+', '-' or None
        :param min_qual: default TopHat: only uniquely mapped list
        :param read_len: only consider reads whose length is in this list
        :return: a dict with starts as key and number of reads as value
        """
        fetch = self.pysam.fetch(chrom, start, stop)
        pos_dict = {}
        for read in fetch:
            if read.mapq >= min_qual:
                if not read_len or read.template_length in read_len:
                    if strand and strand == self.determine_strand(read):
                        pos = read.reference_start
                        if pos not in pos_dict:
                            pos_dict[pos] = 0
                        pos_dict[pos] += 1
        return pos_dict

    def get_coverage(self, chrom, pos, strand='', only_matching=True, min_qual=40):
        """
        get the number of reads at pos
        :param chrom: str chromosome name
        :param pos: int start 0-based
        :param strand: '+', '-' or None
        :param only_matching: reads are only considered in their matching part
        :param min_qual: default TopHat: only uniquely mapped reads
        """
        fetch = self.pysam.fetch(chrom, pos, pos + 1)
        n_reads = 0
        for read in fetch:
            if read.mapq < min_qual:
                continue
            if strand and strand != self.determine_strand(read):
                continue
            if not only_matching or self._type_of_match(read, pos) == 0:
                n_reads += 1
        return n_reads

    def determine_strand(self, read):
        """determines the annotation strand a read would match"""
        if self.reads_orientation == 'mixed':
            return ''
        strand_bool = True
        if read.is_reverse:
            strand_bool = not strand_bool
        if self.reads_orientation == 'reverse':
            strand_bool = not strand_bool
        if read.is_read2:
            strand_bool = not strand_bool
        return '+' if strand_bool else '-'

    @staticmethod
    def del_pos_len(read, max_len=2):
        """
        returns the relative and the absolute position of the first deletion on the read, if any
        :param read: the read (pysam object)
        """
        rel_pos, abs_pos = 0, read.reference_start
        for cigar_token in read.cigartuples:
            if cigar_token[0] == 0:
                rel_pos += cigar_token[1]
                abs_pos += cigar_token[1]
            if cigar_token[0] == 1:
                rel_pos += cigar_token[1]
            if cigar_token[0] == 2 and cigar_token[1] <= max_len:
                return rel_pos, abs_pos, cigar_token[1]
            if cigar_token[0] == 3:
                abs_pos += cigar_token[1]

    @staticmethod
    def _type_of_match(read, pos):
        """
        returns the cigar kind of match between read and reference at pos
        0 -> match, 2 -> del, 3 -> non_match, -1 -> no_overlap
        :param read: the read (pysam object)
        :param pos: ref 0-based position
        """
        tmp_pos = read.reference_start
        if tmp_pos > pos:
            return -1
        cigar_advance = [0,2,3]
        for cigar_token in read.cigartuples:
            if cigar_token[0] in cigar_advance:
                tmp_pos += cigar_token[1]
                if tmp_pos >= pos:
                    return cigar_token[0]
        return -1


class BedGraph:
    """utility functions to parse BedGraph"""

    def __init__(self, path, strand, genome=None, delim='\t', test=False):
        """positions are 0-based"""
        if not test:
            assert os.path.isfile(path)
        self.path = path
        self.strand = strand
        self.delim = delim
        if genome:
            self.chroms_dict = genome.chroms_dict

    def read(self):
        """returns 1-based positions and scores"""
        for linea in open(self.path, 'r'):
            interval = self._parse_line(linea)
            if interval['score'] < 0:
                continue
            for base in range(interval['start'] + 1, interval['stop'] + 1):
                yield base, interval['score']

    def _parse_line(self, linea):
        splat = linea.rstrip('\n').split(self.delim)
        assert len(splat) == 4
        return {'chrom':splat[0], 'start':int(splat[1]), 'stop':int(splat[2]), 'score':float(splat[3])}


class Bed12:
    """utility functions to parse Bed12"""

    def __init__(self, path, delim='', test=False):
        """positions are 0-based"""
        self.path = path
        if not test:
            assert os.path.isfile(path)
            if not delim:
                delim = self._find_delim()
        self.delim = delim

    def _find_delim(self):
        for linea in open(self.path):
            if len(linea.split('\t')) == 12:
                return '\t'
            if len(linea.split(' ')) == 12:
                return ' '
        raise Exception('unknown first line delimiter of bed12 file')

    def read_single_pos(self, bed12_line=None):
        """returns (chrom, pos, strand) for every single pos in each
        interval of the file or of a line, if provided, 0-based"""
        to_iter = iter([bed12_line])
        if not bed12_line:
            to_iter = open(self.path)
        for linea in to_iter:
            parsed = self._parse_line(linea)
            for interval in self._intervals(parsed):
                for pos in range(*interval):
                    yield parsed['chrom'], pos, parsed['strand']

    def _intervals(self, parsed):
        """computes the blocks intervals from a parsed line"""
        for n_block in range(parsed['n_blocks']):
            start = parsed['start'] + parsed['blocks_start'][n_block]
            stop = start+parsed['blocks_size'][n_block]
            yield start, stop
        assert parsed['stop'] == stop

    def _parse_line(self, linea):
        splat = linea.rstrip('\n').split(self.delim)
        assert len(splat) == 12
        return {'chrom':splat[0], 'start':int(splat[1]), 'stop':int(splat[2]), 'name':splat[3], 'score':float(splat[4]),
                'strand':splat[5], 'start_alt':int(splat[6]), 'stop_alt':int(splat[7]), 'n_blocks':int(splat[9]),
                'rgb':[int(x) for x in splat[8].split(',')], 'blocks_size':[int(x) for x in splat[10].split(',') if x],
                'blocks_start':[int(x) for x in splat[11].split(',') if x]}
