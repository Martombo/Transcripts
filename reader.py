import features as ft
import os


class Gtf:
    """utility functions to parse gtf"""

    def __init__(self, gtf_file, gene_ids=None, test=False):
        if not test:
            assert isinstance(gtf_file, str)
            assert os.path.isfile(gtf_file)
        self.gtf_path = gtf_file
        self.gene_ids = gene_ids
        self.genes = {}

    def get_genes(self):
        """get genes dictionary
        :return: genes dict where keys are gene_id
        """
        for linea in open(self.gtf_path).readlines():
            dic = self._parse_exon_line(linea)
            if not dic:
                continue
            gene = self._get_gene(dic)
            if self.gene_ids and dic['attr']['gene_id'] not in self.gene_ids:
                continue
            trans = self._get_transcript(dic, gene)
            exon = ft.Exon(int(dic['attr']['exon_number']), trans, dic['start'], dic['stop'])
        return self.genes

    @staticmethod
    def _get_transcript(dic, gene):
        trans_id = dic['attr']['transcript_id']
        if trans_id in gene.transcripts:
            return gene.transcripts[trans_id]
        else:
            return ft.Transcript(trans_id, gene)

    def _get_gene(self, dic):
        gene_id = dic['attr']['gene_id']
        if gene_id in self.genes:
            return self.genes[gene_id]
        else:
            gene = ft.Gene(gene_id, dic['attr']['gene_name'], dic['chr'], dic['strand'])
            self.genes[gene_id] = gene
            return gene

    def _parse_exon_line(self, linea, feature='.', annot='.'):
        """parses a gtf line, with attributes (see self._parse_attributes())
        :returns {chr, annot, feat, start, ..., attr}"""
        splat = linea.rstrip('\n').split('\t')
        if len(splat) < 8:
            return None
        if splat[2] != 'exon':
            return None
        attr = self._parse_attributes(splat[8])
        return dict(chr=splat[0], annot=splat[1], feat=splat[2], start=int(splat[3]), stop=int(splat[4]),
                    score=splat[5], strand=splat[6], frame=splat[7], attr=attr)

    def _parse_attributes(self, attrs):
        """parses the attributes in a dict
        :returns {'gene_id': 'ENSG000000123', 'gene_version' = '1', ...}
        """
        attrs_splat = attrs.split(';')
        attr_dic = {}
        for attr in attrs_splat:
            if not attr:
                continue
            attr = attr.lstrip(' ')
            attr_key = attr.split(' ')[0]
            attr_item = attr.split('"')[1]
            attr_dic[attr_key] = attr_item
        return attr_dic

