from python.io import hash_fn

class Read(object):

    def __init__(self, seqid, read, min_overlap, direction=0, original_id=None):
        super(Read, self).__init__()
        self.seqid = seqid
        self.read = read
        self.read_length = len(self.read)
        self.direction = direction
        self.matches = []
        self.has_matched = False
        self.min_overlap = min_overlap

        # variables for windowed reads
        self.window_start_idx = None
        self.window_end_idx = None
        self.read_original_length = None
        self.original_id = original_id

        self.r = hash_fn(read)


    def update_matches(self, start, stop):
        self.matches.append([start, stop])
        self.has_matched = True

    def make_windows(self):
        """method for slicing read string into windows

        This is the primary mechanism that facilitates contig extension. If make_windows is not utilized, then
        there will only be alignment within the boundaries of the intiial query and you will get deep coverage instead
        of wide coverage over the query.
        """
        windows = []
        # First create a sliding window read from the left. For example, if we have a string with 5 characters
        # and a minimum overlap of 1 character, then the resulting window indices would be:
        #       |-----| (full Read)
        # [4,5] |---XX|
        # [3,5] |--XXX|
        # [2,5] |-XXXX|
        # [1,5] |XXXXX|
        for i in range(1, self.read_length - self.min_overlap + 1):
            _read = Read("{}_window".format(self.seqid), self.read[i:], self.min_overlap, original_id=self.seqid)
            _read.window_start_idx = i
            _read.window_end_idx = self.read_length
            _read.read_original_length = self.read_length
            windows.append(_read)

        # Next, we want to create a sliding window going in the opposite direction. Using the same example we
        # get the resulting window indices:
        #       |-----| (full Read)
        # [1,2] |XX---|
        # [1,3] |XXX--|
        # [1,4] [XXXX-|
        # [1,5] |XXXXX| (note: this iteration is not run since it was already captured above,
        #                       but is shown here for illustrations

        for i in range(self.min_overlap, self.read_length + 1):
            _read = Read("{}_window".format(self.seqid), self.read[:i], self.min_overlap, original_id=self.seqid)
            _read.window_start_idx = 0
            _read.window_end_idx = i
            _read.read_original_length = self.read_length
            windows.append(_read)


        return windows


class Query(object):

    def __init__(self, q, min_overlap):
        super(Query, self).__init__()
        self.q = q
        self.q_star = q
        self.q_index = list(range(len(q)))
        self.size = self.q_index[-1]
        self.min_overlap = min_overlap

    def make_windows(self, window_size):
        # Slice query into sliding windows of a given size.
        windows = []
        for i in range(self.min_overlap, self.size + window_size - self.min_overlap + 1):
            if i <= self.size:
                if i <= window_size:
                    start = 0
                else:
                    start = i - window_size

                stop = i
            else:
                start = i - window_size
                stop = self.size + 1

            windows.append((start, stop, hash_fn(self.q_star[start:stop])))
        return windows

    def update_q_star(self, left_extension, right_extension):
        """Update q_star and the associated metrics

        left_extension (str): sequence of characters to prepend onto q_star
        right_extension (str): sequence of characters to extend on q_star

        Note: these strings are expected to already be formatted to not contain the characters
        which were shared with the previous q_star. For example, if AABBBB was q_star and CCAA was the read for the
        left extension then we expect the left_extension to be 'CC'. This is expected to be performed using
        contigconstructor.utils.get_extension()
        """
        self.q_star = l_extension + self.q_star + r_extension

        # The index needs to be updated such that the original locations are preserved. Meaning,
        # that the original query always is located at [0, query.size]. Therfore, left extensions will have negative
        # indices and right extensions will have positive indices greater than query.size. We do this here so that
        # when query.make_windows is called the next time, the returned indices will be correct.
        start_idx = self.q_index[0] - len(l_extension)
        stop_idx = self.q_index[1] + len(r_extension)
        self.q_index = list(range(start_idx, stop_idx, 1))
        self.size = self.q_index[-1]



    def hash(self):
        return hash_fn(self.q_star)

