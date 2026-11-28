#!/usr/bin/env python3
"""Export every data table of one or more Terra workspaces to TSV.

Each workspace is given as "<namespace>/<workspace>". Tables are written flat
into --out-dir as "<namespace>.<workspace>.<table>.tsv", in the layout of
Terra's own "Download TSV": an "entity:<table>_id" column followed by the
attributes. List and entity-reference values are written as JSON.
"""
import argparse
import csv
import json
import os
import re
import sys

from terra_notebook_utils import table as tnu_table

# Terra allows only these characters in billing projects, workspaces and table
# names, so "." is a safe separator for the flat file names.
NAME_RE = re.compile(r"^[A-Za-z0-9_-]+$")


def parse_workspace(spec):
    parts = spec.strip().split("/")
    if len(parts) != 2 or not all(NAME_RE.match(p) for p in parts):
        sys.exit(f"Workspace must be '<namespace>/<workspace>' "
                 f"(letters, digits, '_' and '-' only), got: {spec!r}")
    return parts[0], parts[1]


def cell(value):
    if value is None:
        return ""
    if isinstance(value, (list, dict)):
        return json.dumps(value)
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def export_table(namespace, workspace, table_name, out_dir):
    if not NAME_RE.match(table_name):
        sys.exit(f"Unexpected table name {table_name!r} in {namespace}/{workspace}")
    rows = list(tnu_table.list_rows(table_name, workspace=workspace,
                                    workspace_namespace=namespace))
    columns = list(dict.fromkeys(key for row in rows for key in row.attributes))

    path = os.path.join(out_dir, f"{namespace}.{workspace}.{table_name}.tsv")
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow([f"entity:{table_name}_id"] + columns)
        for row in rows:
            writer.writerow([row.name] + [cell(row.attributes.get(c)) for c in columns])
    print(f"{namespace}/{workspace}: {table_name} ({len(rows)} rows, "
          f"{len(columns)} columns) -> {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", action="append", required=True,
                        help="'<namespace>/<workspace>'; repeat for several")
    parser.add_argument("--out-dir", default="terra_tables")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    for spec in dict.fromkeys(args.workspace):
        namespace, workspace = parse_workspace(spec)
        tables = list(tnu_table.list_tables(workspace=workspace,
                                            workspace_namespace=namespace))
        if not tables:
            print(f"{namespace}/{workspace}: no tables")
        for table_name in tables:
            export_table(namespace, workspace, table_name, args.out_dir)


if __name__ == "__main__":
    main()
