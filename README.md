# HomopolyMapper

## Requirements:
python >= 2.7.15
Python libraries: argparse, Bio, json

## Example Use
Compress the reference.fa and generate the map used for lifting between compressed and uncompressed space:

```python build_sparse_compress_map.py -c <path/to/reference.fa> <path/to/output_compressed_reference.fa> <path/to/output_map_file.json> (optional -d <density_of_map>, default=1000)```

If you want to then uncompress the compressed fasta using the map file you just generated:

```python build_sparse_compress_map.py -u <path/to/compressed_reference.fa> <path/to/output_uncompressed_reference.fa> <path/to/map.json>```

If you want to lift a bed file of coordinates from compressed to uncompressed space:

```python lift_seqs.py [-cu , --compressed_to_uncompressed] <path/to/uncompressed_ref.fa> <path/to/map.json> <path/to/compressed_coordinates.bed> <path/to/uncompressed_output.bed>```

If you want to lift a bed file of coordinates from uncompressed to compressed space:

```python lift_seqs.py [-uc , ----uncompressed_to_compressed] <path/to/uncompressed_ref.fa> <path/to/map.json> <path/to/uncompressed_coordinates.bed> <path/to/compressed_output.bed>```
