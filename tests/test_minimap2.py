import unittest
from unittest.mock import patch, mock_open

from scripts.minimap2 import *

class TestCutadapt(unittest.TestCase):
    def test_cutadapt(self):
        # Using mocks to ensure that the cutadapt function is
        # being called correctly without actually running external shell commands
        with patch('subprocess.check_call') as mock_check_call:
            mock_check_call.return_value = 0

            with tempfile.NamedTemporaryFile() as fq1, tempfile.NamedTemporaryFile() as fq2:
                of1_name, of2_name = cutadapt(fq1.name, fq2.name)

            # Check that the mock was called with the correct arguments
            mock_check_call.assert_called_once_with([
                'cutadapt', '-a', 'AGATCGGAAGAGC', '-A', 'AGATCGGAAGAGC',
                '-o', of1_name, '-p', of2_name, '-j', '1', '-m', '10',
                '--quiet', fq1.name, fq2.name
            ])


class TestMakeFilename(unittest.TestCase):
    prefix = "test"
    suffix = "txt"
    def test_make_filename_doesnt_exist(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            filename = make_filename(temp_dir, self.prefix, self.suffix)

            expected_filename = os.path.join(temp_dir, f"{self.prefix}.{self.suffix}")
            self.assertEqual(filename, expected_filename)

    def test_make_filename_without_replace(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            existing_filename = os.path.join(temp_dir, f"{self.prefix}.{self.suffix}")
            with open(existing_filename, 'w') as temp_file:
                temp_file.write("Sample data")

            # Call the function with replace=False
            with patch('sys.exit') as mock_exit, patch('builtins.print') as mock_print:
                filename = make_filename(temp_dir, self.prefix, self.suffix, replace=False)

            expected_filename = os.path.join(temp_dir, f"{self.prefix}.{self.suffix}")
            self.assertEqual(filename, expected_filename)

            # Assert that the function prints the expected message and exits
            mock_print.assert_called_once_with(
                f"Output file {expected_filename} already exists, use --replace to overwrite"
            )
            mock_exit.assert_called_once()

    def test_make_filename_with_replace(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            existing_filename = os.path.join(temp_dir, f"{self.prefix}.{self.suffix}")
            with open(existing_filename, 'w') as temp_file:
                temp_file.write("Sample data")

            # Call the function with replace=True
            with patch('sys.exit') as mock_exit:
                filename = make_filename(temp_dir, self.prefix, self.suffix, replace=True)

            expected_filename = os.path.join(temp_dir, f"{self.prefix}.{self.suffix}")
            self.assertEqual(filename, expected_filename)

            # Assert that the function doesnt exit when replace=True
            mock_exit.assert_not_called()


class TestProcess(unittest.TestCase):
    def setUp(self):
        self.expected_counts = {26272: {'~26273T': {'label': 'aa:E:G10V', 'count': 0.5}},
                                25923: {'~25924T': {'label': 'aa:orf3a:P178S', 'count': 0.5}}}

    def test_process(self):
        counts, coverage = process(fq1="data/unit-tests/minimap2_mutations-coverage_1.fastq",
                         fq2="data/unit-tests/minimap2_mutations-coverage_2.fastq",
                         ref="data/NC_045512.fa", nthread=3, binpath='minimap2', limit=None)
        self.assertEqual(self.expected_counts, counts)

    def test_get_frequencies(self):
        # Zero denominator
        res = [
            {'diff': [('-', 100, 'A')], 'mutations': ['Mutation1']},
            {'diff': [('+', 200, 'C')], 'mutations': ['Mutation2']}
        ]
        coverage = {100: 0, 200: 0}

        result = get_frequencies(res, coverage)
        self.assertEqual(result, {})

        # Non-zero denominator
        res = [
            {'diff': [('-', 100, 'A')], 'mutations': ['Mutation1']},
            {'diff': [('+', 200, 'C')], 'mutations': ['Mutation2']}
        ]
        coverage = {100: 10, 200: 6}

        result = get_frequencies(res, coverage)
        expected_result = {100: {'-101.A': {'count': 0.1, 'label': 'Mutation1'}},
                           200: {'+201.C': {'count': 0.16666666666666666, 'label': 'Mutation2'}}}

        self.assertEqual(result, expected_result)

    def test_write_frequencies(self):
        coverage = {i:0 for i in range(29903)}

        expected_output = 'position,label,mutation,frequency,coverage\n' \
                          '26272,~26273T,aa:E:G10V,0.5,0\n' \
                          '25923,~25924T,aa:orf3a:P178S,0.5,0\n'

        # Mock the open function
        with patch('builtins.open', mock_open()) as mock_file:
            outfile = 'test_output.csv'
            write_frequencies(self.expected_counts, coverage, outfile)

            mock_file.assert_called_once_with(outfile, 'w')

            handle = mock_file()
            calls = handle.write.call_args_list
            written_content = ''.join(call[0][0] for call in calls)

            self.assertEqual(written_content, expected_output)

    def test_write_coverage(self):
        coverage = {
            0: 10,
            1: 6,
            2: 15
        }

        expected_output = 'position,coverage\n' \
                          '0,10\n' \
                          '1,6\n' \
                          '2,15\n'

        # Mock the open function
        with patch('builtins.open', mock_open()) as mock_file:
            outfile = 'test_coverage_output.csv'
            write_coverage(coverage, outfile)

            mock_file.assert_called_once_with(outfile, 'w')

            handle = mock_file()
            calls = handle.write.call_args_list
            written_content = ''.join(call[0][0] for call in calls)

            self.assertEqual(written_content, expected_output)


if __name__ == '__main__':
    unittest.main()