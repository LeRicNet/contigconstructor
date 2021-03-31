import argparse
import numpy as np
import copy

from python.io import read_query_file, logmessage
from python.data import Read, Query
from python.utils import get_extension
from python.veb import VEB

# Parameters
MIN_OVERLAP = 5
QUERY_FILE_PATH = "/Users/ericprince/Documents/CPBS/Courses/CPBS7712/CPBS7712_Module1/Day3Q/data/QUERY.fasta"
TEST_FILE_PATH = "/Users/ericprince/Documents/CPBS/Courses/CPBS7712/CPBS7712_Module1/Day3Q/data/READS.fasta"
SUBSAMPLE = 1000
WINDOW_READS = True
BIDRECTIONAL = True
SAVE_PATH = None
PLOT_CONTIG = True
RANDOM_SEED = 42
VERBOSE = True

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
parser.add_argument("--bidirectional", action='store_true',
                    help="Reads in test_file_path will be read in both directions (L to R and R to L).")
parser.add_argument("--map_mode", type=str, choices=['wide', 'deep'], default="wide")
parser.add_argument("--subsample_n",
                    type=int, required=False, default=0,
                    help="Integer number to sample reads from test_file_path. This can help speed up the algorithm.")
parser.add_argument("--random_seed", type=int, default=1234,
                    help="Seed for random number generator. Only important if subsample_n > 0.")
parser.add_argument("--verbose", action="store_true",
                    help="Print status updates as algorithm is running.")
args, _ = parser.parse_known_args()

print(dir(args))

#
# np.random.seed(RANDOM_SEED)
#
# logmessage("\t ** Loading Query Sequence and Target Sequences", verbose=VERBOSE)
# query = Query(read_query_file(QUERY_FILE_PATH).strip("\n"), MIN_OVERLAP)
# raw_reads = read_query_file(TEST_FILE_PATH, return_all=True)
# reads = [Read(raw_reads[i].strip(">").strip("\n"), raw_reads[i+1].strip("\n"), MIN_OVERLAP, direction=0) for i in range(0, len(raw_reads), 2)]
#
# test_idx = np.random.choice(list(range(len(reads))), SUBSAMPLE)
# reads = [r for i,r in enumerate(reads) if i in test_idx]
#
# rev_reads = copy.deepcopy(reads)
# for r in rev_reads:
#     r.read = r.read[::-1]
#     r.direction = 1
#     r.seqid = "R-{}".format(r.seqid)
#
# windowed_reads = copy.deepcopy(reads)
# rev_windowed_reads = copy.deepcopy(rev_reads)
#
# logmessage("\t ** Making Windowed Reads (F)", verbose=VERBOSE)
# windowed_reads = [r.make_windows() for r in windowed_reads]
# windowed_reads = [item for sublist in windowed_reads for item in sublist]
#
# logmessage("\t ** Making Windowed Reads (R)", verbose=VERBOSE)
# rev_windowed_reads = [r.make_windows() for r in rev_windowed_reads]
# rev_windowed_reads = [item for sublist in rev_windowed_reads for item in sublist]
#
# reads = reads + windowed_reads + rev_windowed_reads + rev_reads
#
# # pickle.dump(reads, open("./reads.p", "wb"))
#
# f_coverage = []
# r_coverage = []
# for i in range(1):
#
#     logmessage("\t ** Creating vEB-scanner", verbose=VERBOSE)
#     # First we build a vEB tree in order to store all of the possible sequences. For each unique read length,
#     # we apply a sliding window across the query (specifically, q_star) in order to capture the edge alignments
#     # as well as any interior matches since read lengths can be less than the query length. For example, a query with
#     # 8 characters, a minimum overlap of 2, and a sliding window of 3 would look like:
#     #
#     #   |--------| (query)
#     #   |XX------| (potential left extension; min overlap)
#     #   |XXX-----| (potential left extension)
#     #   |-XXX----| (potential interior)
#     #   |--XXX---| (potential interior)
#     #   |---XXX--| (potential interior)
#     #   |----XXX-| (potential interior)
#     #   |-----XXX| (potential right extension)
#     #   |------XX| (potential right extension; min overlap)
#     #
#     forward_veb = VEB(u=query.hash() ** 2)
#     unique_read_lengths = np.unique([r.read_length for r in reads])
#     for window_size in unique_read_lengths:
#         windows = query.make_windows(window_size)
#         for window in windows:
#             # window[-1]: hashed string value for window
#             # window[0]: window start location (relative to query.q_index)
#             # window[1]: window end location (relative to query.q_index)
#             forward_veb.insert(window[-1], window[0], window[1])
#
#     # For each read we will check if it is a member in the vEB tree we just created
#     logmessage("\t ** Checking for Target Sequences", verbose=VERBOSE)
#     for read in reads:
#         hits = forward_veb.member(read.r)
#
#         # If there is not a hit in the vEB tree then hits=False. So we only want to look at instances
#         # where hits is not a boolean.
#         if type(hits) != bool:
#
#             # Check if this matched read is a window read (as opposed to the original full read)
#             if 'window' in read.seqid:
#                 # Check for logic that alignment hit is in a feasible location.
#                 # if the read.window_start_idx is 0, then we are only looking for reads which
#                 # extend the query to the right. In this case we want
#                 #
#                 #      hits[0] + read.window_end_idx > query.size
#                 #
#                 # where hits[0] is the starting index for the match within q_star.
#                 if read.window_start_idx == 0 and hits[0] + read.window_end_idx > query.size:
#                     read.update_matches(*hits)
#
#                 # Alternatively, if read.window_end_idx is equal to read.read_original_length then
#                 # we are expecting to read from the left side of the query and are only interested
#                 # in the reads which will extend the query to the left. Accordingly, we want
#                 #
#                 # hits[0] == query.q_index[0]
#                 #
#                 # where hits[0] is the start index of the match and query.q_index[0] is the current right-most index. And,
#                 #
#                 # hits[1] - read.read_length > 0
#                 #
#                 # where hits[1] is the end index for the match within q_star.
#                 if (read.window_end_idx == read.read_original_length and
#                         hits[0] == query.q_index[0] and
#                         hits[1] - read.read_length < 0):
#                     read.update_matches(*hits)
#
#             else:
#                 # If 'window' is not in the read.seqid, then we are considering a full read as a match. These matches
#                 # can only occur within the boundaries of q_star and add depth to the coverage, but no increase in
#                 # length.
#
#                 # TODO: Does this need an additional logical check?
#                 read.update_matches(*hits)
#
#     # Next, we need to organize the reads into matches and mismatches. This is achieved by checking for the has_matched
#     # variable in each Read. We can also separate into forward and reverse matches by checking for 'R' in the seqid.
#     matched_reads = [r for r in reads if r.has_matched]
#     mismatched_reads = [r for r in reads if not r.has_matched]
#
#     forward_matches = [r for r in matched_reads if not r.seqid.startswith('R')]
#     reverse_matches = [r for r in matched_reads if r.seqid.startswith('R')]
#
#     logmessage("\t ** Found {} Forward Orientation Matches, and {} Reverse Orientation Matches.".format(
#         len(forward_matches), len(reverse_matches)), verbose=VERBOSE)
#
#     # Now identify the longest possible extensions on each side using the matches
#     left_extension = {'extension': 0, 'read': None}
#     right_extension = {'extension': 0, 'read': None}
#     interior_matches = []
#     for match in forward_matches + reverse_matches:
#
#         # Extensions have to come from windowed reads in order to get the extension. There can be multiple extensions,
#         # and extensions will often not align with one another so we will just select the largest extension on each side.
#         # These extensions will be matched with any matches found in the interior since that space is fixed for each
#         # iteration. Any reads that were not found in the interior or either of the two largest extensions will be
#         # recycled for another iteration.
#         if 'window' in match.seqid:
#             for _m in match.matches:
#                 # The left extension is identified as having the match start equal to the leftmost query.q_index value and then
#                 # the largest difference between the match location start/stop and the original read_length
#                 if _m[0] == query.q_index[0]:
#                     _extension = match.read_original_length - _m[1]  # _m[1] - _m[0] = _m[1], since _m[0] == 0
#                     if _extension > left_extension['extension']:
#                         _m_best = _m
#                         _r_best = match.read
#                         left_extension['extension'] = _extension
#                         left_extension['read'] = match
#
#                 # The right extension is identified as having a match end equal to the rightmost value of query.q_index and
#                 # similarly having the largest difference between match range and the original read_length
#                 elif _m[1] - 1 == query.q_index[-1]:
#                     _extension = match.read_original_length - (_m[1] - _m[0])
#                     if _extension > right_extension['extension']:
#                         right_extension['extension'] = _extension
#                         right_extension['read'] = match
#
#         else:
#             interior_matches.append(match)
#
#     logmessage("\t ** Extending to the left by {} and to the right by {}".format(
#         left_extension['extension'], right_extension['extension']), verbose=VERBOSE)
#     logmessage("\t ** Found {} interior matches.".format(len(interior_matches)), verbose=VERBOSE)
#
#     # Collect the original reads for each of the extensions
#     print(_r_best)
#     print(_m_best)
#     print(query.q[_m_best[0]:_m_best[1]])
#     tmp = query.q[::-1]
#     print(query.q)
#     print(tmp[_m_best[0]:_m_best[1]])
    # get_extension(left_extension, reads, side="left")

#
#
#     # Separate the reads into matched and mismatched lists
#     f_matched_reads, f_mismatched_reads, r_matched_reads, r_mismatched_reads = [], [], [], []
#     for r in reads:
#         if r.has_matched:
#             if r.direction == 0:
#                 f_matched_reads.append(r)
#             elif r.direction == 1:
#                 r_matched_reads.append(r)
#         else:
#             if r.direction == 0:
#                 f_mismatched_reads.append(r)
#             elif r.direction == 1:
#                 r_mismatched_reads.append(r)
#
#     # matched_reads = [r for r in reads if r.has_matched]
#     print("Found {} matches (F)".format(len(f_matched_reads)))
#     pickle.dump(f_matched_reads, open("./tmp.p", "wb"))
#     matched_windows = [r for r in f_matched_reads if "window" in r.seqid]
#     print("{} windows".format(len(matched_windows)))
#
#     f_max_extension_right, f_max_extension_left = [0, None], [0, None]
#     for matched_read in tqdm(f_matched_reads):
#         _matches = matched_read.matches[-1]
#         if 'window' in matched_read.seqid:
#             original_read = [r for r in reads if r.seqid == matched_read.seqid.split("_")[0]]
#             original_read_size = original_read[0].read_length
#             if _matches[1] - original_read_size < 0:
#                 # print("Extending contig to the left by {} bases".format((int(_matches[1] - original_read_size))))
#                 start = _matches[1] - original_read_size
#                 stop = original_read_size + 1
#                 extension = abs(int(_matches[1] - original_read_size))
#                 if extension > f_max_extension_left[0]:
#                     f_max_extension_left = [extension, matched_read]
#             elif _matches[0] + original_read_size > query.size:
#                 start = _matches[0]
#                 stop = start + original_read_size
#                 extension = abs(int(stop - query.size))
#                 # print("Extending contig to the right by {} bases".format((int(stop - query.size))))
#                 if extension > f_max_extension_right[0]:
#                     f_max_extension_right = [extension, matched_read]
#             else:
#                 start = _matches[0]
#                 stop = _matches[1]
#         else:
#             start = _matches[0]
#             stop = _matches[1]
#
#         f_coverage.extend(list(range(start, stop)))
#
#     # rev_matched_reads = [r for r in rev_reads if r.has_matched]
#     print("Found {} matches (R)".format(len(r_matched_reads)))
#     rev_matched_windows = [r for r in r_matched_reads if "window" in r.seqid]
#     print("{} windows".format(len(rev_matched_windows)))
#
#     r_max_extension_right, r_max_extension_left = [0, None], [0, None]
#     for rev_matched_read in tqdm(r_matched_reads):
#         _matches = rev_matched_read.matches[-1]
#         if 'window' in rev_matched_read.seqid:
#             original_read = [r for r in reads if r.seqid == rev_matched_read.seqid.split("_")[0]]
#             original_read_size = original_read[0].read_length
#             if _matches[1] - original_read_size < 0:
#                 # print("Extending contig to the left by {} bases".format((int(_matches[1] - original_read_size))))
#                 start = _matches[1] - original_read_size
#                 stop = original_read_size + 1
#                 extension = abs(int(_matches[1] - original_read_size))
#                 if extension > r_max_extension_left[0]:
#                     r_max_extension_left = [extension, rev_matched_read]
#             elif _matches[0] + original_read_size > query.size:
#                 start = _matches[0]
#                 stop = start + original_read_size
#                 # print("Extending contig to the right by {} bases".format((int(stop - query.size))))
#                 extension = abs(int(stop - query.size))
#                 if extension > r_max_extension_right[0]:
#                     r_max_extension_right = [extension, rev_matched_read]
#             else:
#                 start = _matches[0]
#                 stop = _matches[1]
#         else:
#             start = _matches[0]
#             stop = _matches[1]
#
#         r_coverage.extend(list(range(start, stop)))
#
#     # TODO: How to determine which contig extension to carry over? Maximum consensus for a given threshold of homology or single maximum?
#
#     right_extension, left_extension = '', ''
#     # Choose the largest extension to create the next query, q_star
#     if f_max_extension_left[1] is not None:
#
#         _id = f_max_extension_left[1].seqid
#         if "R" in _id:
#             _id = _id.split("-")[1]
#         if "window" in _id:
#             _id = _id.split("_")[0]
#         print(_id)
#         # left_extension = [r for r in reads if r.seqid == f_max_extension_left[1].seqid.split("_")[0] and r.direction == 0][0]
#         left_extension = [r for r in reads if r.seqid == _id][0]
#         left_extension = left_extension.read[:-f_max_extension_left[1].read_length]
#
#     if f_max_extension_right[1] is not None:
#         _id = f_max_extension_right[1].seqid
#         if "R" in _id:
#             _id = _id.split("-")[1]
#         if "window" in _id:
#             _id = _id.split("_")[0]
#         right_extension = [r for r in reads if r.seqid == _id][0]
#         right_extension = right_extension.read[f_max_extension_right[1].read_length:]
#
#
#     print("Length L Extension = {}".format(len(left_extension)))
#     print("Length R Extension = {}".format(len(right_extension)))
#     print("Length Q star = {}".format(len(query.q_star)))
#     q_star = left_extension + query.q_star + right_extension
#     print("New length = {}".format(len(q_star)))
#     print(q_star)
#
#     query.update_q_star(left_extension, right_extension)
#     reads = f_mismatched_reads + r_mismatched_reads
#
#     if (f_max_extension_left[1] is None and f_max_extension_right[1] is None) or \
#             (len(reads) == 0):
#         break
#
# # coverage = [list(range(r.matches[-1][0], r.matches[-1][1])) for r in matched]
# #
# # coverage = [item for sublist in coverage for item in sublist]
# #
# fig, (ax1, ax2, ax3) = plt.subplots(nrows=3)
# sns.distplot(f_coverage, color="orange", label="F", ax=ax1)
# sns.distplot(r_coverage, color="blue", label="R", ax=ax2)
# ax3.barh(y=0, height=1, width=query.size, color="gray")
# ax2.set_xlim(ax1.get_xlim())
# ax3.set_xlim(ax1.get_xlim())
#
# with open("./qstar.txt", "w") as f:
#     f.write(query.q_star)
#
#
# plt.show()


