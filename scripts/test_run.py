print("Python being used:")
import sys
print(sys.executable)

print("\nTesting core dependencies...\n")

import numpy as np
print("numpy OK")

from ase import Atoms
print("ase OK")

import calorine
print("calorine OK")

import polaritonic_deep_md
print("polaritonic_deep_md OK:", polaritonic_deep_md.__file__)

print("\nALL TESTS PASSED")