![alt text](examples/images/readfish_logo.jpg "ReadFish Logo")

Configuration Files
===

Specification for the TOML files used to drive readfish experiments.

---

Table of Contents
===
 - [TOML files](#toml-files)
 - [Config sections](#config-sections)
   - [Guppy connection](#guppy-connection)
   - [Conditions](#conditions)
 
 
TOML files
===
ReadFish experiments are configured using TOML files, which are minimal and 
easy-to-read markup files. Each line is either a key-value pair or a 'table' 
heading. See more in the [TOML specification](https://github.com/toml-lang/toml).


Config sections
===

Guppy connection
---
The `caller_settings` table specifies the basecalling parameters used by guppy.  

The `config_name` parameter must a valid guppy configuration excluding the file 
extension; these can be found in the `data` folder of the your guppy installation 
directory (`/opt/ont/guppy/data/*.cfg`).  

### Remote basecalling

```toml
[caller_settings]
config_name = "dna_r9.4.1_450bps_fast"
host = "REMOTE_SERVER_IP_ADDRESS"
port = "REMOTE_GUPPY_SERVER_PORT"
``` 

### Local basecalling

```toml
[caller_settings]
config_name = "dna_r10.4.1_e8.2_400bps_fast"
host = "ipc:///tmp/.guppy"
port = 5555
```

### Barcoding

Although Readfish enables per-barcode targeting, GNASTy has not be formally tested using barcodes. See the main [Readfish](https://github.com/LooseLab/readfish) repository if you wish to try this.

Conditions
---

If using a graph index for graph alignment, specify this here:
```toml
[conditions]
reference = "/absolute/path/to/reference.gfa"
```

### Conditions sub-tables

Each conditions sub-table must contain all of the following keys, these are the same between barcoded and non-barcoded toml files.

|          Key |       Type      | Values | Description |
|-------------:|:---------------:|:------:|:------------|
|name|string|N/A|The name given to this condition|
|control|bool|N/A|Is this a control condition. If `true` all reads will be ignored in this region|
|min_chunks|int|N/A|The minimum number of read chunks to evaluate|
|max_chunks|int|N/A|The maximum number of read chunks to evaluate|
|targets|string or array|N/A|The genomic targets to accept or reject; see [types](#target-types) and [formats](#target-formats)|
|single_on|string|[unblock, stop_receiving, proceed]|The action to take when a read has a single on-target mapping|
|multi_on|string|[unblock, stop_receiving, proceed]|The action to take when a read has multiple on-target mappings|
|single_off|string|[unblock, stop_receiving, proceed]|The action to take when a read has a single off-target mapping|
|multi_off|string|[unblock, stop_receiving, proceed]|The action to take when a read has multiple off-target mappings|
|no_seq|string|[unblock, stop_receiving, proceed]|The action to take when a read does not basecall|
|no_map|string|[unblock, stop_receiving, proceed]|The action to take when a read does not map to your reference|

The physical layout of each flowcell constrains how many experimental conditions 
can be used; the number of sub-tables in the `conditions` section determines how 
the flowcell is divided. 

The maximum number of conditions for MinION and PromethION flowcells is given in 
the table below. The number of conditions must be a factor of the number for the
selected combination.

<table>
  <tr>
    <td rowspan="2"></td>
    <th colspan="2">Axis</th>
  </tr>
  <tr>
    <th>0</th>
    <th>1</th>
  </tr>
  <tr>
    <th>MinION</th>
    <td>16</td>
    <td>32</td>
  </tr>
  <tr>
    <th>PromethION</th>
    <td>25</td>
    <td>120</td>
  </tr>
</table>

#### When using graph alignment

For graph alignment, the following parameters are ignored, as they are overwritten by the `--len_cutoff` and `--align_threshold` command line arguments, or are not available using a graph.
- `targets`
- `multi_on`
- `multi_off`
- `no_seq`

You should also set `min_chunks = 0` and `max_chunks = 0`, as read length is controlled by `--len_cutoff`.

See the [example](https://github.com/samhorsfield96/readfish/blob/graph_alignment_bifrost/examples/human_chr_selection.toml) TOML file for running with GNASTy.
