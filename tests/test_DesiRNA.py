import unittest
import tempfile
from DesiRNA import initialize_simulation
from DesiRNA import InputFile, Nucleotide, read_input

class TestInitializeSimulation(unittest.TestCase):

    def test_initialize_with_standard_design(self):
        # Create an InputFile instance with test data
        test_input = InputFile(
            name="Design",
            sec_struct="((((((.((((((((....))))).)).).))))))",
            seq_restr="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
        )
        nucleotides = initialize_simulation(test_input)
        
        self.assertIsInstance(nucleotides, list)  # Check it's a list
        self.assertTrue(all(isinstance(n, Nucleotide) for n in nucleotides))  # Check all items are Nucleotide instances
        self.assertEqual(len(nucleotides), 36)  # Check the list has the correct length

        # Check properties of a specific Nucleotide object
        # Example: Check the eighth Nucleotide in the list
        nucleotide_8 = nucleotides[8]
        self.assertEqual(nucleotide_8.number, 8)
        self.assertEqual(sorted(nucleotide_8.letters), sorted(['A', 'C', 'G', 'U']))
        self.assertEqual(nucleotide_8.pairs_with, 26)
        self.assertEqual(sorted(nucleotide_8.pair_letters), sorted(['U', 'C', 'G', 'A']))
        self.assertEqual(sorted(nucleotide_8.letters_allowed), sorted(['U', 'A', 'G', 'C']))


class TestReadInput(unittest.TestCase):

    def create_temp_input_file(self, content):
        temp_file = tempfile.NamedTemporaryFile(mode='w', delete=False)
        temp_file.write(content)
        temp_file.close()
        return temp_file.name

    def test_standard_input(self):
        # Prepare mock file content
        mock_file_content = ">name\nDesign\n>sec_struct\n((((((.((((((((....))))).)).).))))))\n>seq_restr\nNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n"
        temp_file_path = self.create_temp_input_file(mock_file_content)

        # Call read_input and test results
        input_file = read_input(temp_file_path)
        self.assertEqual(input_file.name, "Desin")
        self.assertEqual(input_file.sec_struct, "((((((.((((((((....))))).)).).))))))")
        self.assertEqual(input_file.seq_restr, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")


if __name__ == '__main__':
    unittest.main()
