version 1.0

# Builds a Map from each of Keys[j] to [Values[0][j], Values[1][j], ...] --
# i.e. transposes Values (each sub-array a "column" aligned with Keys) into
# one entry per key gathering that key's value across all columns. E.g.
# Keys = ["A", "B", "C"], Values = [["1","2","3"], ["4","5","6"]] produces
# {"A": ["1","4"], "B": ["2","5"], "C": ["3","6"]}.
task BuildFileMap {
  input {
    Array[String] Keys
    Array[Array[String]] Values
    String ImageTag = "latest"
    Int MemoryGB = 4
    Int? DiskGB
  }

  File keys_json = write_json(Keys)
  File values_json = write_json(Values)

  command <<<
    set -euo pipefail

    python3 <<CODE
import json

with open("~{keys_json}") as f:
    keys = json.load(f)
with open("~{values_json}") as f:
    columns = json.load(f)

result = {key: [col[j] for col in columns] for j, key in enumerate(keys)}
json.dump(result, open("manifest.json", "w"))
CODE
  >>>

  runtime {
    docker: "ayenkin1871/mei-lr-association-python_general:" + ImageTag
    memory: MemoryGB + " GB"
    cpu: 2
    disks: "local-disk " + select_first([DiskGB, 20]) + " SSD"
    preemptible: 3
    maxRetries: 2
  }

  output {
    Map[String, Array[File]] FileMap = read_json("manifest.json")
  }
}
