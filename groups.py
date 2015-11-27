import reader as rd


class TranscriptsGroup:

    def __init__(self, gtf=None):
        self.transcripts = {}
        if gtf:
            assert isinstance(gtf, rd.Gtf)
            self.gtf = gtf

    def ordered_transcripts_tss(self):
        if not self.gtf.genes:
            self.gtf.get_genome()
        for gene in self.gtf.genes.values():
            for transcript in gene.transcripts.values():
                self._add_tss(transcript)

    def get_antisense_transcripts(self, trans_list=None):
        antisenses = []
        for trans_chrom in self.transcripts.values():
            for trans_i in range(len(trans_chrom)):
                trans = trans_chrom[trans_i]
                if trans_list and trans.id not in trans_list:
                    continue
                antisenses.append([trans, self._get_closest_antisense(trans, trans_i, trans_chrom)])
        return antisenses

    def _get_closest_antisense(self, trans, trans_i, trans_chrom):
        strand = trans.strand
        as_downstream = self._as_in_range(range(trans_i + 1, len(trans_chrom)), trans_chrom, strand, trans.gene.id)
        as_upstream = self._as_in_range(range(trans_i - 1, -1, -1), trans_chrom, strand, trans.gene.id)
        if not as_upstream:
            return as_downstream
        if not as_downstream:
            return as_upstream
        if abs(trans.tss - as_downstream.tss) < abs(trans.tss - as_upstream.tss):
            return as_downstream
        return as_upstream

    @staticmethod
    def _as_in_range(trans_range, trans_chrom, strand, gene_id):
        for trans_j in trans_range:
            if trans_chrom[trans_j] != strand and trans_chrom[trans_j].gene.id != gene_id:
                return trans_chrom[trans_j]

    def _add_tss(self, trans):
        if trans.chromosome.name not in self.transcripts:
            self.transcripts[trans.chromosome.name] = [trans]
        else:
            chrom_list = self.transcripts[trans.chromosome.name]
            insert_pos = len(chrom_list)
            for upstream_trans in chrom_list[::-1]:
                if upstream_trans.tss < trans.tss:
                    chrom_list[insert_pos:insert_pos] = [trans]
                    break
                insert_pos -= 1
            else:
                chrom_list[0:0] = [trans]
