"""Validate existing exported components, without importing their claimed proofs.

Run from the repository root: python test/check_artifact_schemas.py
Requires jsonschema and Python 3.11+; does not modify artifacts.
"""
import json
from pathlib import Path
import tomllib

from jsonschema import Draft202012Validator
from referencing import Registry, Resource

root = Path(__file__).resolve().parents[1]
registry = Registry()
schema_urls = {}
for version in (1, 2):
    path = root / "schemas" / f"qubo-component-v{version}.schema.json"
    schema = json.loads(path.read_text(encoding="utf-8"))
    Draft202012Validator.check_schema(schema)
    registry = registry.with_resource(path.as_uri(), Resource.from_contents(schema))
    schema_urls[f"qubo-component/{version}"] = path.as_uri()

count = 0
for directory in (root / "perf/q2/results/patterns", root / "perf/q3/results", root / "perf/q4/results", root / "perf/q5/results", root / "perf/q6/results", root / "perf/q7/results", root / "perf/q9/results"):
    for path in directory.rglob("*.toml"):
        with path.open("rb") as source:
            data = tomllib.load(source)
        version = data.get("schema_version")
        if version in schema_urls:
            Draft202012Validator({"$ref": schema_urls[version]}, registry=registry).validate(data)
            count += 1
print(f"Validated {count} component artifacts against their versioned schemas")

# Q8 campaign records wrap their component. Runtime benchmark caches and bound
# envelopes are deliberately not counted as additional scientific constructions.
wrapped_count = 0
for path in (root / "perf/q8/results/records").glob("*.toml"):
    with path.open("rb") as source:
        record = tomllib.load(source)
    data = record.get("component")
    if data is None:
        continue
    version = data["schema_version"]
    assert version in schema_urls, (path, version)
    Draft202012Validator({"$ref": schema_urls[version]}, registry=registry).validate(data)
    wrapped_count += 1
print(f"Validated {wrapped_count} Q8 wrapped component artifacts (not necessarily certified PASS)")
