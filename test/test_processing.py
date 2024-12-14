import unittest
import snarl_analyser

class TestFileProcessor(unittest.TestCase):
    
    def test_count_words(self):
        # Temporary file for testing
        test_content = "hello world hello"
        with open("test.txt", "w") as f:
            print(test_content)

if __name__ == "__main__":
    unittest.main()

