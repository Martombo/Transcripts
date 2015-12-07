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

    def get_genome(self):
        """returns a Genome object, from an Ensembl gtf file"""
        with open(self.gtf_path) as f:
            for linea in f.readlines():
                line_dic = self._parse_line(linea)
                gene = self._get_gene(line_dic)
                if not line_dic or 'transcript_id' not in line_dic or \
                        (self.gene_biotype and line_dic['gene_biotype'] not in self.gene_biotype):
                    continue
                trans = self._get_transcript(line_dic, gene)
                if line_dic['type'] == 'start_codon':
                    trans.cds_start = min(line_dic['start'], line_dic['stop'])
                elif line_dic['type'] == 'stop_codon':
                    trans.cds_stop = max(line_dic['start'], line_dic['stop'])
                elif line_dic['type'] == 'exon':
                    exon = ft.Exon(int(line_dic['exon_number']), trans, line_dic['start'], line_dic['stop'])
                elif line_dic['type'] == 'CDS':
                    trans.add_cds(line_dic['start'], line_dic['stop'])
        return self.genome

    @staticmethod
    def _get_transcript(attr, gene):
        """
        returns a Transcript from the Gene, if it's already present.
        otherwise a new Transcript
        """
        trans_id = attr['transcript_id']
        if trans_id in gene.trans_dict:
            return gene.trans_dict[trans_id]
        return ft.Transcript(trans_id, gene)

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

    def get_coverage(self, chrom, start, stop=None, strand='', min_qual=40):
        """
        get the number of reads in region
        reads is counted even if only 1 base overlaps region
        :param chrom: str chromosome name
        :param start: int start 0-based
        :param stop: int stop 0-based
        :param strand: '+', '-' or None
        :param min_qual: default TopHat: only uniquely mapped reads
        """
        if not stop:
            stop = start + 1
        fetch = self.pysam.fetch(chrom, start, stop)
        n_reads = 0
        for read in fetch:
            if read.mapq >= min_qual:
                if strand and strand == self.determine_strand(read):
                    n_reads += 1
        return n_reads

    def determine_strand(self, read):
        """determines the annotation strand a read would match"""
        if self.reads_orientation == 'mixed':
            return 'NA'
        strand_bool = True
        if read.is_reverse:
            strand_bool = not strand_bool
        if self.reads_orientation == 'reverse':
            strand_bool = not strand_bool
        if read.is_read2:
            strand_bool = not strand_bool
        return '+' if strand_bool else '-'


class BedGraph:
    """utility functions to parse BedGraph"""

    def __init__(self, path, strand, genome=None, test=False):
        """positions are 0-based"""
        if not test:
            assert os.path.isfile(path)
        self.path = path
        self.strand = strand
        if genome:
            self.chroms_dict = genome.chroms_dict

    def iter(self):
        """returns 1-based positions and scores"""
        for linea in open(self.path, 'r'):
            interval = self._parse_line(linea)
            if interval['score'] < 0:
                continue
            for base in range(interval['start'] + 1, interval['stop']):
                yield base, interval['score']


    @staticmethod
    def _parse_line(linea, delim='\t'):
        splat = linea.rstrip('\n').split(delim)
        assert len(splat) == 4
        return {'chrom': splat[0], 'start': int(splat[1]), 'stop': int(splat[2]), 'score':float(splat[3])}
