import unittest as ut
import unittest.mock as um
import features as ft
import reader as rd
import pysam
import os


class TestFeatureGeneral(ut.TestCase):
    def test_fix_order(self):
        sites = ft.fix_order(10, 20, '+')
        self.assertLess(sites[0], sites[1])

    def test_fix_order_opp(self):
        sites = ft.fix_order(20, 10, '+')
        self.assertLess(sites[0], sites[1])

    def test_fix_order_rev(self):
        sites = ft.fix_order(10, 20, '-')
        self.assertGreater(sites[0], sites[1])

    def test_fix_order_rev_opp(self):
        sites = ft.fix_order(20, 10, '-')
        self.assertGreater(sites[0], sites[1])

    def test_move_pos(self):
        pos = ft.move_pos(100, +10, '+')
        self.assertEqual(110, pos)

    def test_move_pos2(self):
        pos = ft.move_pos(100, +10, '-')
        self.assertEqual(90, pos)

    def test_move_pos3(self):
        pos = ft.move_pos(100, -10, '-')
        self.assertEqual(110, pos)

    def test_before_strand(self):
        self.assertTrue(ft.before_strand(100, 200, '+'))

    def test_before_strand_rev(self):
        self.assertTrue(ft.before_strand(200, 100, '-'))

    def test_before_strand_false(self):
        self.assertFalse(ft.before_strand(200, 100, '+'))

    def test_motif_to_sequences_simple(self):
        motif = ft.Motif('ACGT')
        sequences = motif.possible_sequences
        self.assertEquals(1, len(sequences))
        self.assertIn('ACGT', sequences)

    def test_motif_to_sequences_1base_2choices(self):
        motif = ft.Motif('ARGT')
        sequences = motif.possible_sequences
        self.assertEquals(2, len(sequences))
        self.assertIn('AAGT', sequences)
        self.assertIn('AGGT', sequences)

    def test_motif_to_sequences_2bases_2choices(self):
        motif = ft.Motif('ARGTK')
        sequences = motif.possible_sequences
        self.assertEquals(4, len(sequences))
        self.assertIn('AAGTG', sequences)
        self.assertIn('AAGTT', sequences)
        self.assertIn('AGGTG', sequences)
        self.assertIn('AGGTT', sequences)

    def test_motif_to_sequences_36(self):
        motif = ft.Motif('RACATTDGCCHTGAGY')
        sequences = motif.possible_sequences
        self.assertEquals(36, len(sequences))


class TestExons(ut.TestCase):
    chr1 = ft.Chromosome('chr1')
    gene1 = ft.Gene('g1', chr1, 'gene1', '+')
    gene1_rev = ft.Gene('g1', chr1, 'gene1', '-')

    def test_distance_from_start(self):
        trans1 = ft.Transcript('t1', self.gene1)
        exon1 = ft.Exon('e1', 1, trans1, 100, 200)
        self.assertEqual(10, exon1.relative_position(110))

    def test_distance_from_start_rev(self):
        trans1 = ft.Transcript('t1', self.gene1_rev)
        exon1 = ft.Exon('e1', 1, trans1, 100, 200)
        self.assertEqual(70, exon1.relative_position(130))

    def test_includes(self):
        trans1 = ft.Transcript('t1', self.gene1_rev)
        exon1 = ft.Exon('e1', 1, trans1, 100, 200)
        self.assertTrue(exon1.includes(exon1.start))
        self.assertTrue(exon1.includes(exon1.genomic_stop))

    def test_sequence(self):
        genom = ft.Genome(fasta_path='/Users/martin/notDropbox/utils/gNome/Homo_sapiens.GRCh38.dna.primary_assembly.fa')
        chr1 = ft.Chromosome('1', genome=genom)
        gene1 = ft.Gene('g1', chr1, 'gene1', '+')
        trans1 = ft.Transcript('t1', gene1)
        start = 11869
        stop = 12227
        exon1 = ft.Exon('e1', 1, trans1, start, stop)
        seq = exon1.get_sequence()
        self.assertEqual(stop - start + 1, len(seq))
        self.assertEqual('GTTAACTTGCCGTCAGCCTTTTCTTTGACCTCTTCTTTCTGTTCATGTGTATTTGCTGTCTCTTAGCCCAGACTTCCCGTGTCCTTTCC' +
                         'ACCGGGCCTTTGAGAGGTCACAGGGTCTTGATGCTGTGGTCTTCATCTGCAGGTGTCTGACTTCCAGCAACTGCTGGCCTGTGCCAGGG' +
                         'TGCAAGCTGAGCACTGGAGTGGAGTTTTCCTGTGGAGAGGAGCCATGCCTAGAGTGGGATGGGCCATTGTTCATCTTCTGGCCCCTGTT' +
                         'GTCTGCATGTAACTTAATACCACAACCAGGCATAGGGGAAAGATTGGAGGAAAGATGAGTGAGAGCATCAACTTCTCTCACAACCTAGG' +
                         'CCA', seq)


class TestGene(ut.TestCase):
    chr1 = ft.Chromosome('chr1')
    gene1 = ft.Gene('g1', chr1, 'gene1', '+')
    trans1 = ft.Transcript('t1', gene1)
    trans2 = ft.Transcript('t2', gene1)
    exon1 = ft.Exon('e1', 1, trans1, 100, 200)
    exon2 = ft.Exon('e2', 2, trans1, 300, 400)
    exon3 = ft.Exon('e3', 1, trans2, 500, 700)
    exon4 = ft.Exon('e4', 2, trans2, 800, 900)

    def test_transcripts(self):
        n = 0
        for trans in self.gene1.transcripts:
            n += 1
        self.assertEqual(2, n)

    def test_included1(self):
        self.assertTrue(self.gene1.includes(150))

    def test_included2(self):
        self.assertTrue(self.gene1.includes(850))

    def test_not_included1(self):
        self.assertFalse(self.gene1.includes(250))

    def test_not_included2(self):
        self.assertFalse(self.gene1.includes(1000))


class TestTranscripts(ut.TestCase):
    chr1 = ft.Chromosome('chr1')
    gene1 = ft.Gene('g1', chr1, 'gene1', '+')
    gene2 = ft.Gene('g2', chr1, 'gene2', '-')
    trans1 = ft.Transcript('t1', gene1)
    ft.Exon('e1', 1, trans1, 1, 300)
    ft.Exon('e2', 2, trans1, 501, 1000)
    trans1.add_cds(101, 800)
    trans2 = ft.Transcript('t2', gene2)
    ft.Exon('e1', 1, trans2, 1300, 801)
    ft.Exon('e2', 2, trans2, 600, 501)
    ft.Exon('e3', 3, trans2, 300, 101)
    trans2.add_cds(1200, 200)
    m = um.Mock()

    def test_no_splice_sites(self):
        self.assertEquals(0, len(self.trans1.splice_sites))

    def test_2exons(self):
        trans1 = ft.Transcript('t1', self.gene1)
        ft.Exon('e1', 1, trans1, 100, 200)
        ft.Exon('e2', 2, trans1, 300, 400)
        self.assertEquals(2, len(trans1.exon_dict))

    def test_2trans(self):
        ft.Transcript('t2', self.gene1)
        self.assertEquals(2, len(list(self.gene1.trans_dict.values())))

    def test_add_cds(self):
        self.assertEqual(101, self.trans1.cds_start)
        self.assertEqual(800, self.trans1.cds_stop)

    def test_add_cds_rev(self):
        trans = ft.Transcript('t', self.gene2)
        trans.add_cds(100, 200)
        self.assertEqual(200, trans.cds_start)
        self.assertEqual(100, trans.cds_stop)

    def test_add_2cds(self):
        trans = ft.Transcript('t', self.gene1)
        trans.add_cds(100, 200)
        trans.add_cds(500, 600)
        self.assertEqual(100, trans.cds_start)
        self.assertEqual(600, trans.cds_stop)

    def test_add_2cds_rev(self):
        trans = ft.Transcript('t', self.gene2)
        trans.add_cds(500, 600)
        trans.add_cds(100, 200)
        self.assertEqual(600, trans.cds_start)
        self.assertEqual(100, trans.cds_stop)

    def test_abs_pos_dw_1exon(self):
        trans = ft.Transcript('t', self.gene1)
        ft.Exon('e1', 1, trans, 100, 200)
        pos = trans.abs_pos_downstream(150, 10)
        self.assertEqual(160, pos)

    def test_abs_pos_dw_2exons(self):
        trans = ft.Transcript('t', self.gene1)
        ft.Exon('e1', 1, trans, 100, 200)
        ft.Exon('e2', 2, trans, 300, 400)
        pos = trans.abs_pos_downstream(150, 100)
        self.assertEqual(349, pos)

    def test_abs_pos_dw_3exons(self):
        trans = ft.Transcript('t', self.gene1)
        ft.Exon('e1', 1, trans, 100, 200)
        ft.Exon('e2', 2, trans, 300, 400)
        ft.Exon('e3', 3, trans, 800, 900)
        pos = trans.abs_pos_downstream(150, 200)
        self.assertEqual(848, pos)

    def test_abs_pos_dw_1exon_rev(self):
        trans = ft.Transcript('t', self.gene2)
        ft.Exon('e1', 1, trans, 100, 200)
        pos = trans.abs_pos_downstream(150, 10)
        self.assertEqual(140, pos)

    def test_abs_pos_dw_2exons_rev(self):
        trans = ft.Transcript('t', self.gene2)
        ft.Exon('e1', 1, trans, 300, 400)
        ft.Exon('e2', 2, trans, 100, 200)
        pos = trans.abs_pos_downstream(350, 100)
        self.assertEqual(151, pos)

    def test_abs_pos_dw_3exons_rev(self):
        trans = ft.Transcript('t', self.gene2)
        ft.Exon('e1', 1, trans, 600, 700)
        ft.Exon('e2', 2, trans, 300, 400)
        ft.Exon('e3', 3, trans, 100, 200)
        pos = trans.abs_pos_downstream(650, 200)
        self.assertEqual(152, pos)

    def test_abs_pos_up_1exon(self):
        trans = ft.Transcript('t', self.gene1)
        ft.Exon('e1', 1, trans, 100, 200)
        pos = trans.abs_pos_upstream(150, 10)
        self.assertEqual(140, pos)

    def test_abs_pos_up_2exons(self):
        trans = ft.Transcript('t', self.gene1)
        ft.Exon('e1', 1, trans, 100, 200)
        ft.Exon('e2', 2, trans, 300, 400)
        pos = trans.abs_pos_upstream(340, 100)
        self.assertEqual(141, pos)

    def test_abs_pos_up_3exons(self):
        trans = ft.Transcript('t', self.gene1)
        ft.Exon('e1', 1, trans, 100, 200)
        ft.Exon('e2', 2, trans, 300, 400)
        ft.Exon('e3', 3, trans, 800, 923)
        pos = trans.abs_pos_upstream(820, 210)
        self.assertEqual(112, pos)

    def test_abs_pos_up_1exon_rev(self):
        trans = ft.Transcript('t', self.gene2)
        ft.Exon('e1', 1, trans, 100, 200)
        pos = trans.abs_pos_upstream(150, 10)
        self.assertEqual(160, pos)

    def test_abs_pos_up_2exons_rev(self):
        trans = ft.Transcript('t', self.gene2)
        ft.Exon('e1', 1, trans, 400, 500)
        ft.Exon('e2', 2, trans, 100, 200)
        pos = trans.abs_pos_upstream(120, 110)
        self.assertEqual(429, pos)

    def test_intervals_1exon(self):
        trans = ft.Transcript('t', self.gene1)
        ft.Exon('e1', 1, trans, 100, 200)
        inter = list(trans.intervals(140, 160))[0]
        self.assertEqual(140, inter[0])
        self.assertEqual(160, inter[1])

    def test_intervals_2exon(self):
        trans = ft.Transcript('t', self.gene1)
        ft.Exon('e1', 1, trans, 100, 200)
        ft.Exon('e2', 2, trans, 300, 400)
        inter = list(trans.intervals(140, 360))
        self.assertEqual(140, inter[0][0])
        self.assertEqual(200, inter[0][1])
        self.assertEqual(300, inter[1][0])
        self.assertEqual(360, inter[1][1])

    def test_intervals_2exon_rev(self):
        trans = ft.Transcript('t', self.gene2)
        ft.Exon('e1', 1, trans, 300, 400)
        ft.Exon('e2', 2, trans, 100, 200)
        inter = list(trans.intervals(360, 140))
        self.assertEqual(360, inter[0][0])
        self.assertEqual(300, inter[0][1])
        self.assertEqual(200, inter[1][0])
        self.assertEqual(140, inter[1][1])

    def test_cds_relative_starts_simple(self):
        self.m.get_read_starts = lambda *args: {199:1}
        rel_starts = self.trans1.get_cds_relative_starts(self.m)
        self.assertEquals(100, list(rel_starts.keys())[0])

    def test_cds_relative_starts_multiple(self):
        self.m.get_read_starts = lambda *args: {199:1, 249:1}
        rel_starts = self.trans1.get_cds_relative_starts(self.m)
        keys = list(rel_starts.keys())
        self.assertEquals(2, len(keys))
        self.assertEquals(100, keys[0])
        self.assertEquals(150, keys[1])

    def test_cds_relative_starts_first(self):
        self.m.get_read_starts = lambda *args: {100:1}
        rel_starts = self.trans1.get_cds_relative_starts(self.m)
        self.assertEquals(1, list(rel_starts.keys())[0])

    def test_cds_relative_starts_before(self):
        self.m.get_read_starts = lambda *args: {99:1}
        rel_starts = self.trans1.get_cds_relative_starts(self.m)
        self.assertEquals(0, len(list(rel_starts.keys())))

    def test_cds_relative_starts_exon(self):
        self.m.get_read_starts = lambda *args: {599:1}
        rel_starts = self.trans1.get_cds_relative_starts(self.m)
        self.assertEquals(300, list(rel_starts.keys())[0])

    def test_cds_relative_starts_rev_simple(self):
        self.m.get_read_starts = lambda *args: {1100:1}
        rel_starts = self.trans2.get_cds_relative_starts(self.m)
        self.assertEquals(100, list(rel_starts.keys())[0])

    def test_cds_relative_starts_rev_first(self):
        self.m.get_read_starts = lambda *args: {1199:1}
        rel_starts = self.trans2.get_cds_relative_starts(self.m)
        self.assertEquals(1, list(rel_starts.keys())[0])

    def test_cds_relative_starts_rev_exon(self):
        self.m.get_read_starts = lambda *args: {599:1}
        rel_starts = self.trans2.get_cds_relative_starts(self.m)
        self.assertEquals(401, list(rel_starts.keys())[0])

    def test_cds_relative_starts_rev_exons(self):
        self.m.get_read_starts = lambda *args: {299:1}
        rel_starts = self.trans2.get_cds_relative_starts(self.m)
        self.assertEquals(501, list(rel_starts.keys())[0])


class TestChromosome(ut.TestCase):
    chr1 = ft.Chromosome('chr1')
    gene1 = ft.Gene('g1', chr1)
    trans1 = ft.Transcript('t1', gene1)
    trans2 = ft.Transcript('t2', gene1)
    gene2 = ft.Gene('g2', chr1)
    trans3 = ft.Transcript('t3', gene2)

    def test_transcripts(self):
        n = 0
        for k in self.chr1.transcripts:
            n += 1
        self.assertEqual(3, n)


class TestSequence(ut.TestCase):
    def test_is_start(self):
        seq = ft.Sequence('AAAAAATGCCCCCC')
        self.assertFalse(seq.is_start(0))
        self.assertFalse(seq.is_start(11))
        self.assertTrue(seq.is_start(5))

    def test_is_stop(self):
        seq = ft.Sequence('AAAATAGCCGCCCC')
        self.assertFalse(seq.is_stop(0))
        self.assertFalse(seq.is_stop(10))
        self.assertTrue(seq.is_stop(4))

    def test_next_stop(self):
        seq = ft.Sequence('ACGACGTGACCGCCCC')
        self.assertEqual(6, seq.next_stop(0))
        self.assertIsNone(seq.next_stop(1))
        self.assertIsNone(seq.next_stop(7))

    def test_get_orfs(self):
        seq = ft.Sequence('AATGCCCTGACCC')
        orfs = seq.get_orfs()
        self.assertEqual(1, len(orfs))
        self.assertEqual(1, orfs[0][0])
        self.assertEqual(7, orfs[0][1])

    def test_get_2orfs(self):
        seq = ft.Sequence('AATGCCCTGACCCGATGTAA')
        orfs = seq.get_orfs()
        self.assertEqual(2, len(orfs))
        self.assertEqual(1, orfs[0][0])
        self.assertEqual(7, orfs[0][1])
        self.assertEqual(14, orfs[1][0])
        self.assertEqual(17, orfs[1][1])

    def test_rnaplfold(self):
        seq = ft.Sequence('aaaaaaaaagggggggggttttttttt')
        folds = seq.run_RNAplfold()
        self.assertGreater(folds, 0)


class TestReaderGtf(ut.TestCase):
    attr = 'gene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_id "ENST00000249857"; exon_number "1";\n'
    gtf_line = '\t'.join(['1', 'havana', 'exon', '11869', '14409', '.', '+', '.', attr])

    def test_parse_line(self):
        gtf = rd.Gtf('', test=True)
        line_dic = gtf._parse_line(self.gtf_line)
        self.assertIsNotNone(line_dic)
        gene = gtf._get_gene(line_dic)
        self.assertIsNotNone(gene)

    def test_it(self):
        gtf_file = '/Users/martin/notDropbox/utils/genes/Homo_sapiens.GRCh38.81_head.gtf'
        if os.path.isfile(gtf_file):
            gtf = rd.Gtf('/Users/martin/notDropbox/utils/genes/Homo_sapiens.GRCh38.81_head.gtf')
            genome = gtf.get_genome()
            self.assertIsNotNone(genome)


class TestReaderBam(ut.TestCase):
    read = pysam.AlignedSegment()
    parser_bam = rd.Bam(test=True, reads_orientation='forward')
    read.reference_start = 100
    read.cigartuples = [[0,10],[3,10],[0,10],[1,10],[0,10]]

    def test_det_strand(self):
        self.read.is_reverse = False
        self.read.is_read2 = False
        strand = self.parser_bam.determine_strand(self.read)
        self.assertEquals('+', strand)

    def test_det_strand_rev(self):
        self.read.is_reverse = True
        self.read.is_read2 = False
        self.assertNotEqual('reverse', self.parser_bam.reads_orientation)
        strand = self.parser_bam.determine_strand(self.read)
        self.assertEquals('-', strand)

    def test_det_strand_rev_mate(self):
        self.read.is_reverse = True
        self.read.is_read2 = True
        strand = self.parser_bam.determine_strand(self.read)
        self.assertEquals('+', strand)

    def test_det_strand_opp(self):
        self.read.is_reverse = False
        self.read.is_read2 = False
        self.parser_bam.reads_orientation = 'reverse'
        strand = self.parser_bam.determine_strand(self.read)
        self.parser_bam.reads_orientation = 'forward'
        self.assertEquals('-', strand)

    def test_type_of_match(self):
        type_of_match = self.parser_bam._type_of_match(self.read, 105)
        self.assertEqual(0, type_of_match)

    def test_type_of_match1(self):
        type_of_match = self.parser_bam._type_of_match(self.read, 115)
        self.assertEqual(3, type_of_match)

    def test_type_of_match2(self):
        type_of_match = self.parser_bam._type_of_match(self.read, 125)
        self.assertEqual(0, type_of_match)

    def test_type_of_match3(self):
        type_of_match = self.parser_bam._type_of_match(self.read, 135)
        self.assertEqual(0, type_of_match)


class TestReaderBed12(ut.TestCase):
    parser_bed12 = rd.Bed12('', '', test=True)
    bed12_line = '5	100	350	gene1	0	+	100	350	0	3	10,30,50,	0,100,200,\n'
    bed12_line_space = '5 100 350 gene1 0 + 100 350 0 3 10,30,50 0,100,200,\n'
    m = um.Mock()

    def test_find_delim(self):
        self.m.return_value.__iter__ = lambda x: iter([self.bed12_line])
        with um.patch('reader.open', self.m, create=True):
            self.assertEquals('\t', self.parser_bed12._find_delim())

    def test_find_delim_space(self):
        self.m.return_value.__iter__ = lambda x: iter([self.bed12_line_space])
        with um.patch('reader.open', self.m, create=True):
            self.assertEquals(' ', self.parser_bed12._find_delim())

    def test_parse_line(self):
        self.parser_bed12.delim = '\t'
        parsed = self.parser_bed12._parse_line(self.bed12_line)
        self.assertEquals(12, len(parsed))
        self.assertEquals('gene1', parsed['name'])
        self.assertEquals('+', parsed['strand'])
        self.assertEquals([10,30,50], parsed['blocks_size'])

    def test_intervals(self):
        self.parser_bed12.delim = '\t'
        parsed = self.parser_bed12._parse_line(self.bed12_line)
        intervals = list(self.parser_bed12._intervals(parsed))
        self.assertEquals(3, len(intervals))
        self.assertEquals((200,230), intervals[1])

    def test_read_single_pos_line(self):
        self.parser_bed12.delim = '\t'
        single_poss = list(self.parser_bed12.read_single_pos(self.bed12_line))
        self.assertEquals(('5', 100, '+'), single_poss[0])
        self.assertEquals(90, len(single_poss))

    def test_read_single_pos(self):
        self.parser_bed12.delim = '\t'
        self.m.return_value.__iter__ = lambda x: iter([self.bed12_line])
        with um.patch('reader.open', self.m, create=True):
            single_poss = list(self.parser_bed12.read_single_pos())
            self.assertEquals(('5', 100, '+'), single_poss[0])
            self.assertEquals(90, len(single_poss))
