import features as ft
import os


class Gtf:
    """utility functions to parse gtf"""

    def __init__(self, gtf_file, gene_ids=None, test=False):
        """
        :param gtf_file: path to gtf file
        :param gene_ids: all genes whose id is not in this list are discarded
        :param test: for testing
        """
        if not test:
            assert isinstance(gtf_file, str)
            assert os.path.isfile(gtf_file)
        self.gtf_path = gtf_file
        self.gene_ids = gene_ids
        self.genes = {}

    def get_genes(self):
        """get genes dictionary
        :returns genes dict where keys are gene_id
        """
        for linea in open(self.gtf_path).readlines():
            line_dic = self._parse_line(linea)
            if not line_dic:
                continue
            gene = self._get_gene(line_dic)
            if self.gene_ids and line_dic['gene_id'] not in self.gene_ids:
                continue
            if line_dic['type'] == 'stop_codon':
                trans = self._get_transcript(line_dic, gene)
                trans.cds_stop = max(line_dic['start'], line_dic['stop'])
            elif line_dic['type'] == 'exon':
                trans = self._get_transcript(line_dic, gene)
                exon = ft.Exon(int(line_dic['exon_number']), trans, line_dic['start'], line_dic['stop'])
        return self.genes

    @staticmethod
    def _get_transcript(attr, gene):
        trans_id = attr['transcript_id']
        if trans_id in gene.transcripts:
            return gene.transcripts[trans_id]
        else:
            return ft.Transcript(trans_id, gene)

    def _get_gene(self, dic):
        gene_id = dic['gene_id']
        if gene_id in self.genes:
            return self.genes[gene_id]
        else:
            gene = ft.Gene(gene_id, dic['gene_name'], dic['chr'], dic['strand'])
            self.genes[gene_id] = gene
            return gene

    def _parse_line(self, string):
        """parses a gtf line, with attributes (see self._parse_attributes())
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
