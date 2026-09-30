# Tests

This directory contains automated tests for the package.

Run the complete test suite from the repository root:

```bash
pytest -v
```

Tests are being developed from low-level components toward complete workflows:

```text
imports
  ↓
geometry
  ↓
score parsing and filtering
  ↓
glycan extraction
  ↓
library construction
  ↓
placement
  ↓
orientation analysis
  ↓
end-to-end regression
```

Small intentional fixtures belong under `tests/data/`. Large modeling datasets and manual regression calculations should remain under the Git-ignored `work/` directory.
