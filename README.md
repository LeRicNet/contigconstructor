# ContigConstructor

## Motivation

Contigs are genomic sequences
We want to see what the surrounding milieu is with respect to a known sequence of interest.

## Methodology


### Dependencies

ContigConstructor was specifically built under Python 3.6.8 and uses NumPy 1.17.2. While this can not be guaranteed, we expect that ContigConstructor can be run on Python 3+.

## Easy Install

Assuming dependencies have been met, ContigConstructor can utilized with the following simple commands:

```
cd <DESTNATION-DIR>
git clone https://github.com/princeew/contigconstructor
cd ./contigconstructor
```

## Usage

The main script for ContigConstructor is `./build_contig.py`. It can be executed from the command line using the following specifications:

```
build_contig.py [-h] --query_file_path QUERY_FILE_PATH --test_file_path
                       TEST_FILE_PATH --save_path SAVE_PATH [--bidirectional]
                       [--map_mode {wide,deep}] [--subsample_n SUBSAMPLE_N]
                       [--random_seed RANDOM_SEED] [--verbose]

optional arguments:
  -h, --help            show this help message and exit
  --query_file_path QUERY_FILE_PATH
                        Path to FASTA formatted file containing single query
                        sequence. This sequence is the starting seed for
                        contig growth.
  --test_file_path TEST_FILE_PATH
                        Path to FASTA formatted file containing sequencing
                        reads. These reads are what will be mapped to the
                        query sequence. NOTE: Reads are expected to be
                        preprocessed for trimming adapter sequences.
  --save_path SAVE_PATH
                        Path to directory where output files will be saved.
  --bidirectional       Reads in test_file_path will be read in both
                        directions (L to R and R to L).
  --map_mode {wide,deep}
  --subsample_n SUBSAMPLE_N
                        Integer number to sample reads from test_file_path.
                        This can help speed up the algorithm.
  --random_seed RANDOM_SEED
                        Seed for random number generator. Only important if
                        subsample_n > 0.
  --verbose             Print status updates as algorithm is running.
  
```

Note that `--map_mode=wide` reflects that both the query sequence 

<p align="center">
  <img src="./examples/imgs/map_mode.svg">
</p>

### Input File Format

Both `QUERY_FILE_PATH` and `TEST_FILE_PATH` are required to be in FASTA format. For example the query file might look like like the following:

```
>INITIAL_QUERY
GGGATCGGCCATTGAACAAGATGGATTGCACGCAGGTT
```

And the file provided to `TEST_FILE_PATH` might look like:

```
>2S43D:03629:08794
TTCAGGCTCTGGCATGCATTAGAAATGTGGCTTG
>2S43D:08938:01257
GGGTGGTCCCCCTCCTTTACTTGTAACGTTGTCCTAAG
```

### Outputs


## Example

The following command will execute a subsampled run on example data:

```
python build_contig.py \
--query_file_path=./examples/data/QUERY.fasta  \
--test_file_path=./examples/data/READS.fasta \
--save_path=./examples/contig/ \
--bidirectional \
--map_mode=wide \
--subsample_n=1000 \
--verbose
```

Example outputs can be found under `./examples/contig/`.

## Limitations

### Low-Complexity Regions

### Sequence Read Quality

### Phase 1 Limitations

The current implementation (Phase 1) only performs one pass over the reads provided to TEST_FILE_PATH,

## Troubleshooting

Please report any issues, comments, or questions to Eric Prince via email Eric.Prince@CUAnschutz.edu, or [file an issue](https://github.com/princeew/contigconstructor/issues).
