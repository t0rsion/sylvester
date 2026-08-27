import os
import unittest

import driver


class PortableCommandTests(unittest.TestCase):
    def test_harness_path_becomes_relative(self):
        argument = os.path.join(driver.HERE, "inputs", "cyclic-4.sing")

        self.assertEqual(
            driver.portable_command(["Singular", argument]),
            ["Singular", os.path.join("inputs", "cyclic-4.sing")],
        )

    def test_external_path_keeps_only_the_file_name(self):
        argument = os.path.abspath(os.path.join(driver.HERE, "..", "private", "tool"))

        self.assertEqual(driver.portable_command([argument]), ["tool"])


if __name__ == "__main__":
    unittest.main()
