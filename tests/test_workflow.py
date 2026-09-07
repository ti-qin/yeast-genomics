"""Data-free regression tests; fake executables do not validate biology."""
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'workflow' / 'lib')]
from rna_io import extract_counts, validate_manifest
from find_read import find_read
from sequence_utils import getNonNBED, Fasta

class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='yeast test ')
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        self.bin = self.work / 'bin'
        self.bin.mkdir()
        self.env = dict(os.environ, PATH=str(self.bin) + os.pathsep + os.environ['PATH'],
                        PYTHONDONTWRITEBYTECODE='1')

    def fake(self, name, body):
        path = self.bin / name
        path.write_text('#!/usr/bin/env python3\n' + body)
        path.chmod(0o755)

    def run_script(self, name, *args, success=True):
        result = subprocess.run(['bash', str(ROOT / 'workflow' / name), *map(str,args)], cwd=self.work,
                                env=self.env, text=True, capture_output=True, timeout=30)
        if success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def test_shell_syntax_and_help(self):
        for script in sorted((ROOT / 'workflow').glob('*.sh')):
            with self.subTest(script=script.name):
                subprocess.run(['bash','-n',str(script)],check=True)
                self.run_script(script.name,'--help')

    def test_missing_option_and_zero_jobs(self):
        result=self.run_script('03_assembly_nextdenovo.sh','-r',success=False)
        self.assertIn('requires a value',result.stderr)
        result=self.run_script('03_assembly_nextdenovo.sh','-j','0',success=False)
        self.assertIn('positive integer',result.stderr)

    def test_single_n_and_coordinate_convention(self):
        for sequence, length in [('ANAA',3),('N',0),('NNAAAN',3),('AAAA',4),('',0)]:
            self.assertEqual(len(getNonNBED('contig',sequence)),length)
        bed=getNonNBED('contig','ANAA')
        self.assertEqual([(x.start,x.end) for x in bed.coordinates[0]],[(1,2),(3,5)])

    def test_mito_requires_two_evidence_classes_and_respects_length(self):
        # Coverage-only, homology-only and protein-only must remain nuclear.
        data = "\n".join([
            "coverage\t10000\t80\t20\t1\t0\t0",
            "dna\t10000\t20\t20\t0\t1\t0",
            "protein\t10000\t20\t20\t0\t0\t3",
            "coverage_dna\t200000\t80\t20\t1\t1\t0",
            "dna_protein\t10000\t20\t20\t0\t1\t3",
            "below_protein_threshold\t10000\t20\t20\t0\t1\t2",
            "too_long\t200001\t80\t20\t1\t5\t3",
        ]) + "\n"
        result = subprocess.run(['awk','-v','m=200000','-f',str(ROOT/'workflow/lib/classify_mito.awk')],
                                input=data, text=True, capture_output=True, check=True)
        flags = {row.split('\t')[0]: int(row.split('\t')[-1]) for row in result.stdout.splitlines()}
        self.assertEqual([name for name,flag in flags.items() if flag], ['coverage_dna','dna_protein'])

    def test_empty_fasta(self):
        path=self.work/'empty.fa';path.touch()
        self.assertEqual(len(Fasta(str(path))),0)

    def test_read_matching_is_exact_and_unambiguous(self):
        reads=self.work/'reads';reads.mkdir()
        for name in ['G10_1.fq.gz','G1_1.fq.gz','G1_2.fq.gz']:(reads/name).touch()
        self.assertEqual(find_read(reads,'G1',r'[_.-]1\.fq\.gz$').name,'G1_1.fq.gz')
        (reads/'G1_lane2_1.fq.gz').touch()
        with self.assertRaises(ValueError):find_read(reads,'G1',r'[_.-]1\.fq\.gz$')
        with self.assertRaises(ValueError):find_read(reads,'missing',r'[_.-]1\.fq\.gz$')

    def test_count_header_keeps_first_gene(self):
        raw=self.work/'counts.txt';out=self.work/'matrix.tsv'
        raw.write_text('# Program:featureCounts\nGeneid\tChr\tStart\tEnd\tStrand\tLength\t/path/sample.sorted.bam\ngene1\tchr1\t1\t9\t+\t9\t4\ngene2\tchr1\t20\t30\t+\t11\t7\n')
        extract_counts(raw,out)
        self.assertEqual(out.read_text(),'Geneid\tsample\ngene1\t4\ngene2\t7\n')

    def test_manifest_rejects_conflicting_references_and_duplicate_samples(self):
        files=[self.work/name for name in ['x_R1.fq.gz','x_R2.fq.gz','genome.fa','genes.gtf','other.fa']]
        for path in files:path.write_text('fixture\n')
        row=['yeast',*map(str,files[:4])]
        manifest=self.work/'runs.tsv';manifest.write_text('\t'.join(row)+'\n')
        validate_manifest(manifest)
        manifest.write_text(('\t'.join(row)+'\n')*2)
        with self.assertRaisesRegex(ValueError,'duplicate'):validate_manifest(manifest)
        other=row.copy();other[3]=str(files[4]);manifest.write_text('\t'.join(row)+'\n'+'\t'.join(other)+'\n')
        with self.assertRaisesRegex(ValueError,'inconsistent'):validate_manifest(manifest)

    def test_nextdenovo_isolated_output_and_worker_failure(self):
        reads=self.work/'reads';reads.mkdir();(reads/'G1.fastq.gz').touch()
        self.fake('nextDenovo', '''import configparser, pathlib, sys
cfg=configparser.ConfigParser();cfg.read(sys.argv[1])
out=pathlib.Path(cfg['General']['workdir'])/'03.ctg_graph'
out.mkdir(parents=True)
(out/'nd.asm.fasta').write_text('>ctg\\nAAAA\\n')
''')
        self.run_script('03_assembly_nextdenovo.sh','-r',reads,'-w',self.work/'output','-j','1')
        self.assertTrue((self.work/'output/G1.nextdenovo.raw.fasta').is_file())
        self.fake('nextDenovo','import sys\nsys.exit(17)\n')
        self.run_script('03_assembly_nextdenovo.sh','-r',reads,'-w',self.work/'failed','-j','1',success=False)

    def test_nextpolish_output_and_failed_worker(self):
        inputs=self.work/'assemblies';inputs.mkdir();(inputs/'G1.recon3.fasta').write_text('>c\nAA\n')
        reads=self.work/'reads';reads.mkdir()
        for mate in (1,2):(reads/f'G1_{mate}.fq.gz').write_text('fixture\n')
        self.fake('nextPolish', "from pathlib import Path\nPath('genome.nextpolish.fasta').write_text('>c\\nAA\\n')\n")
        self.run_script('08_polish_nextpolish.sh','-a',inputs,'-r',reads,'-o',self.work/'out','-j','1')
        self.assertTrue((self.work/'out/G1.nextpolish.fasta').is_file())
        self.fake('nextPolish','import sys\nsys.exit(11)\n')
        self.run_script('08_polish_nextpolish.sh','-a',inputs,'-r',reads,'-o',self.work/'failed','-j','1',success=False)

    def test_necat_does_not_continue_after_correction_failure(self):
        reads=self.work/'reads';reads.mkdir();(reads/'G1.fastq.gz').touch()
        self.fake('necat.pl', '''import sys
assert sys.argv[1] == 'correct', 'continued after failure'
sys.exit(19)
''')
        result=self.run_script('02_assembly_necat.sh','-r',reads,'-o',self.work/'out','-j','1',success=False)
        self.assertNotIn('All done',result.stdout)

    def test_compleasm_summary_and_failure(self):
        inputs=self.work/'assemblies';inputs.mkdir();(inputs/'G1.fasta').write_text('>c\nAA\n')
        library=self.work/'library';library.mkdir()
        self.fake('compleasm', '''import pathlib, sys
assert sys.argv[1]=='run'
out=pathlib.Path(sys.argv[sys.argv.index('-o')+1]);out.mkdir(parents=True)
(out/'summary.txt').write_text('S: 90.0%\\nD: 5.0%\\nF: 2.0%\\nM: 3.0%\\n')
''')
        self.run_script('10_assess_compleasm.sh','-i',inputs,'-o',self.work/'out','-L',library,'-j','1')
        table=(self.work/'out/BuscoScores.tsv').read_text()
        self.assertIn('G1\t95\t90.0\t5.0\t2.0\t3.0',table)
        self.fake('compleasm','import sys\nsys.exit(1)\n')
        self.run_script('10_assess_compleasm.sh','-i',inputs,'-o',self.work/'failed','-L',library,'-j','1',success=False)

    def test_blast_failure_does_not_produce_filtered_assembly(self):
        draft=self.work/'draft.fa';draft.write_text('>c\nAAAA\n')
        self.fake('blastn','import sys\nsys.exit(1)\n')
        result=subprocess.run([sys.executable,str(ROOT/'workflow/lib/remove_redundant_contigs.py'),'-d',str(draft),'-o',str(self.work/'filtered')],cwd=self.work,env=self.env,capture_output=True)
        self.assertNotEqual(result.returncode,0)
        self.assertFalse((self.work/'filtered.NR.fasta').exists())

    def test_redundancy_threshold_80_and_single_n(self):
        draft=self.work/'draft.fa';draft.write_text('>a\nAAAANAAAAA\n>b\nAAAAAAAAAA\n')
        self.fake('blastn', '''import pathlib,sys
pathlib.Path(sys.argv[sys.argv.index('-out')+1]).write_text('a\\tb\\t100\\t9\\t0\\t0\\t1\\t9\\t1\\t9\\t0\\t50\\n')
''')
        result=subprocess.run([sys.executable,str(ROOT/'workflow/lib/remove_redundant_contigs.py'),'-d',str(draft),'-o',str(self.work/'filtered')],cwd=self.work,env=self.env,capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(Fasta(str(self.work/'filtered.NR.fasta')).getID(),['b'])

    def test_qv_dry_run_handles_spaces_without_execution(self):
        inputs=self.work/'assemblies';inputs.mkdir();(inputs/'G1.nextpolish.fasta').write_text('>a\nAA\n')
        reads=self.work/'reads';reads.mkdir()
        for mate in (1,2):(reads/f'G1_{mate}.fq.gz').touch()
        result=self.run_script('09_assess_merqury.sh','-a',inputs,'-r',reads,'-o',self.work/'out','--dry-run')
        self.assertIn('[DRY-RUN]',result.stdout)
        self.assertFalse((self.work/'out/Collapsed.qv').exists())

if __name__ == '__main__':
    unittest.main()
