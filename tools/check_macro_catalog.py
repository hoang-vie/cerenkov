#!/usr/bin/env python3
"""Static policy check for user-facing macro organization.

The simulator deliberately keeps all .mac files flat because CMake exposes
macros/*.mac in the build directory.  This checker ensures each macro has a
recognized intent category and prevents public examples from drifting back to
legacy aliases or omitting the minimum run/source/readout/output structure.
"""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
MACRO_DIR = ROOT / "macros"

SPECIAL = {
    "init_vis.mac": "visualization",
    "vis.mac": "visualization",
}

PREFIXES = {
    "example_": "example",
    "validation_": "validation",
    "study_": "study",
    "r1_shell_": "legacy",
}

LEGACY_IN_PUBLIC_EXAMPLES = {
    "/detector/coreMaterial",
    "/detector/sipmSize",
    "/detector/sipmPitch",
    "/detector/setAcrylicThickness",
    "/detector/useAcrylicShell",
    "/detector/setSiPMOffsetY",
}

REQUIRED_EXAMPLE_TOKENS = {
    "/detector/coreProfile",
    "/detector/sensorMode",
    "/output/fileName",
    "/run/initialize",
    "/generator/mode",
    "/run/beamOn",
}


def classify(name: str):
    if name in SPECIAL:
        return SPECIAL[name]
    for prefix, category in PREFIXES.items():
        if name.startswith(prefix):
            return category
    return None


def active_lines(path: Path):
    lines = []
    for raw in path.read_text(errors="replace").splitlines():
        line = raw.split("#", 1)[0].strip()
        if line:
            lines.append(line)
    return lines


def main() -> int:
    errors = []
    counts = {"example": 0, "validation": 0, "study": 0, "legacy": 0, "visualization": 0}
    macros = sorted(MACRO_DIR.glob("*.mac"))

    for path in macros:
        category = classify(path.name)
        if category is None:
            errors.append(f"uncategorized macro: macros/{path.name}")
            continue
        counts[category] += 1

        if category == "example":
            lines = active_lines(path)
            text = "\n".join(lines)
            for legacy in sorted(LEGACY_IN_PUBLIC_EXAMPLES):
                if any(line.startswith(legacy) for line in lines):
                    errors.append(f"public example uses legacy command {legacy}: macros/{path.name}")
            missing = sorted(tok for tok in REQUIRED_EXAMPLE_TOKENS if tok not in text)
            if missing:
                errors.append(
                    f"public example missing required structure {', '.join(missing)}: macros/{path.name}"
                )

    readme = MACRO_DIR / "README.md"
    if not readme.exists():
        errors.append("missing macros/README.md")
    else:
        readme_text = readme.read_text(errors="replace")
        for path in macros:
            if classify(path.name) == "example" and path.name not in readme_text:
                errors.append(f"public example not catalogued in macros/README.md: {path.name}")

    if errors:
        print("Macro catalog checks: FAIL")
        for err in errors:
            print(f"  - {err}")
        return 1

    print(
        "Macro catalog checks: OK "
        f"({counts['example']} examples; {counts['validation']} validation; "
        f"{counts['study']} studies; {counts['legacy']} legacy; "
        f"{counts['visualization']} visualization)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
