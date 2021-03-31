import argparse
import copy
import numpy as np
import os
import pandas as pd

from python.io import read_query_file, logmessage
from python.data import Read, Query
from python.veb import VEB

# Load User Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--query_file_path",
                    type=str, required=True,
                    help="Path to FASTA formatted file containing single query sequence. "
                         "This sequence is the starting seed for contig growth.")
parser.add_argument("--test_file_path",
                    type=str, required=True,
                    help="Path to FASTA formatted file containing sequencing reads. "
                         "These reads are what will be mapped to the query sequence. "
                         "NOTE: Reads are expected to be preprocessed for trimming adapter sequences.")
parser.add_argument("--save_path",
                    type=str, required=True,
                    help="Path to directory where output files will be saved.")
parser.add_argument("--min_overlap", type=int, required=True, default=5,
                    help="Minimum number of bases that must overlap for a windowed read to match.")
parser.add_argument("--bidirectional", action='store_true',
                    help="Reads in test_file_path will be read in both directions (L to R and R to L).")
parser.add_argument("--map_mode", type=str, choices=['wide', 'deep'], default="wide",
                    help="Whether to extend the contig by creating windowed reads, or only map full matches within the query sequence.")
parser.add_argument("--subsample_n",
                    type=int, required=False, default=0,
                    help="Integer number to sample reads from test_file_path. This can help speed up the algorithm.")
parser.add_argument("--random_seed", type=int, default=1234,
                    help="Seed for random number generator. Only important if subsample_n > 0.")
parser.add_argument("--verbose", action="store_true",
                    help="Print status updates as algorithm is running.")
args, _ = parser.parse_known_args()

# Set Parameters
QUERY_FILE_PATH = args.query_file_path
TEST_FILE_PATH = args.test_file_path
MIN_OVERLAP = args.min_overlap
SUBSAMPLE = args.subsample_n
WINDOW_READS = True if args.map_mode == "wide" else False
BIDRECTIONAL = args.bidirectional
SAVE_PATH = args.save_path
RANDOM_SEED = args.random_seed
VERBOSE = args.verbose
# PLOT_CONTIG = Not Implemented

np.random.seed(RANDOM_SEED)

if not os.path.isdir(SAVE_PATH):
    os.makedirs(SAVE_PATH)

# Load Data
logmessage("\t ** Loading Query Sequence and Target Sequences", verbose=VERBOSE)
query = Query(read_query_file(QUERY_FILE_PATH).strip("\n"), MIN_OVERLAP)
raw_reads = read_query_file(TEST_FILE_PATH, return_all=True)
reads = [Read(raw_reads[i].strip(">").strip("\n"), raw_reads[i+1].strip("\n"), MIN_OVERLAP, direction=0) for i in range(0, len(raw_reads), 2)]

test_idx = np.random.choice(list(range(len(reads))), SUBSAMPLE)
reads = [r for i,r in enumerate(reads) if i in test_idx]

if BIDRECTIONAL:
    rev_reads = copy.deepcopy(reads)
    for r in rev_reads:
        r.read = r.read[::-1]
        r.direction = 1
        r.seqid = "R-{}".format(r.seqid)


windowed_reads = copy.deepcopy(reads)
rev_windowed_reads = copy.deepcopy(rev_reads)

logmessage("\t ** Making Windowed Reads (F)", verbose=VERBOSE)
windowed_reads = [r.make_windows() for r in windowed_reads]
windowed_reads = [item for sublist in windowed_reads for item in sublist]

logmessage("\t ** Making Windowed Reads (R)", verbose=VERBOSE)
rev_windowed_reads = [r.make_windows() for r in rev_windowed_reads]
rev_windowed_reads = [item for sublist in rev_windowed_reads for item in sublist]

if BIDRECTIONAL and WINDOW_READS:
    reads = reads + windowed_reads + rev_windowed_reads + rev_reads
#
# # pickle.dump(reads, open("./reads.p", "wb"))
#
f_coverage = []
r_coverage = []
# for i in range(1):
#
logmessage("\t ** Creating vEB-scanner", verbose=VERBOSE)
# First we build a vEB tree in order to store all of the possible sequences. For each unique read length,
# we apply a sliding window across the query (specifically, q_star) in order to capture the edge alignments
# as well as any interior matches since read lengths can be less than the query length. For example, a query with
# 8 characters, a minimum overlap of 2, and a sliding window of 3 would look like:
#
#   |--------| (query)
#   |XX------| (potential left extension; min overlap)
#   |XXX-----| (potential left extension)
#   |-XXX----| (potential interior)
#   |--XXX---| (potential interior)
#   |---XXX--| (potential interior)
#   |----XXX-| (potential interior)
#   |-----XXX| (potential right extension)
#   |------XX| (potential right extension; min overlap)
#
forward_veb = VEB(u=query.hash() ** 2)
unique_read_lengths = np.unique([r.read_length for r in reads])
for window_size in unique_read_lengths:
    windows = query.make_windows(window_size)
    for window in windows:
        # window[-1]: hashed string value for window
        # window[0]: window start location (relative to query.q_index)
        # window[1]: window end location (relative to query.q_index)
        forward_veb.insert(window[-1], window[0], window[1])

# For each read we will check if it is a member in the vEB tree we just created
logmessage("\t ** Checking for Target Sequences", verbose=VERBOSE)
for read in reads:
    hits = forward_veb.member(read.r)

    # If there is not a hit in the vEB tree then hits=False. So we only want to look at instances
    # where hits is not a boolean.
    if type(hits) != bool:

        # Check if this matched read is a window read (as opposed to the original full read)
        if 'window' in read.seqid:
            # Check for logic that alignment hit is in a feasible location.
            # if the read.window_start_idx is 0, then we are only looking for reads which
            # extend the query to the right. In this case we want
            #
            #      hits[0] + read.window_end_idx > query.size
            #
            # where hits[0] is the starting index for the match within q_star.
            if read.window_start_idx == 0 and hits[0] + read.window_end_idx > query.size:
                read.update_matches(*hits)

            # Alternatively, if read.window_end_idx is equal to read.read_original_length then
            # we are expecting to read from the left side of the query and are only interested
            # in the reads which will extend the query to the left. Accordingly, we want
            #
            # hits[0] == query.q_index[0]
            #
            # where hits[0] is the start index of the match and query.q_index[0] is the current right-most index. And,
            #
            # hits[1] - read.read_length > 0
            #
            # where hits[1] is the end index for the match within q_star.
            if (read.window_end_idx == read.read_original_length and
                    hits[0] == query.q_index[0] and
                    hits[1] - read.read_length < 0):
                read.update_matches(*hits)

        else:
            # If 'window' is not in the read.seqid, then we are considering a full read as a match. These matches
            # can only occur within the boundaries of q_star and add depth to the coverage, but no increase in
            # length.

            # TODO: Does this need an additional logical check?
            read.update_matches(*hits)

# Next, we need to organize the reads into matches and mismatches. This is achieved by checking for the has_matched
# variable in each Read. We can also separate into forward and reverse matches by checking for 'R' in the seqid.
matched_reads = [r for r in reads if r.has_matched]
mismatched_reads = [r for r in reads if not r.has_matched]

forward_matches = [r for r in matched_reads if not r.seqid.startswith('R')]
reverse_matches = [r for r in matched_reads if r.seqid.startswith('R')]

logmessage("\t ** Found {} Forward Orientation Matches, and {} Reverse Orientation Matches.".format(
    len(forward_matches), len(reverse_matches)), verbose=VERBOSE)

# Now identify the longest possible extensions on each side using the matches
left_extension = {'extension': 0, 'read': None}
right_extension = {'extension': 0, 'read': None}
interior_matches = []
for match in forward_matches + reverse_matches:

    # Extensions have to come from windowed reads in order to get the extension. There can be multiple extensions,
    # and extensions will often not align with one another so we will just select the largest extension on each side.
    # These extensions will be matched with any matches found in the interior since that space is fixed for each
    # iteration. Any reads that were not found in the interior or either of the two largest extensions will be
    # recycled for another iteration.
    if 'window' in match.seqid:
        for _m in match.matches:
            # The left extension is identified as having the match start equal to the leftmost query.q_index value and then
            # the largest difference between the match location start/stop and the original read_length
            if _m[0] == query.q_index[0]:
                _extension = match.read_original_length - _m[1]  # _m[1] - _m[0] = _m[1], since _m[0] == 0
                if _extension > left_extension['extension']:
                    _m_best = _m
                    _r_best = match.read
                    left_extension['extension'] = _extension
                    left_extension['read'] = match

            # The right extension is identified as having a match end equal to the rightmost value of query.q_index and
            # similarly having the largest difference between match range and the original read_length
            elif _m[1] - 1 == query.q_index[-1]:
                _extension = match.read_original_length - (_m[1] - _m[0])
                if _extension > right_extension['extension']:
                    right_extension['extension'] = _extension
                    right_extension['read'] = match

    else:
        interior_matches.append(match)

logmessage("\t ** Extending to the left by {} and to the right by {}".format(
    left_extension['extension'], right_extension['extension']), verbose=VERBOSE)
logmessage("\t ** Found {} interior matches.".format(len(interior_matches)), verbose=VERBOSE)

# Phase 1 Output
# The desired output is to provide seqids with index information
all_matches = [left_extension['read']] + interior_matches + [right_extension['read']]
results = {
    'sseqid': [],
    'qseqid': [],
    'sstart': [],
    'ssend': [],
    'qstart': [],
    'qend': [],
    'direction': []
}
for _m in all_matches:
    results['sseqid'].append(_m.original_id)
    results['qseqid'].append('contig:1')  # Current Phase 1 implementation will only result in one contig
    results['sstart'].append(_m.matches[-1][0])
    results['ssend'].append(_m.matches[-1][1])
    results['qstart'].append(_m.window_start_idx)
    results['qend'].append(_m.window_end_idx)
    results['direction'].append("F" if not _m.seqid.startswith("R") else "R")

pd.DataFrame(results).to_csv(os.path.join(SAVE_PATH, "output.aln"))