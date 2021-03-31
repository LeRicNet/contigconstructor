"""General utility functions for contigconstructor"""


def get_extension(extension_dict, reads, side):
    """Isolate extended sequence for q_star update

    extension_dict (dict): extension information including the read and size of the extension
    reads (list): a list of contigconstructor.data.Read objects
    side (str): either 'left' or 'right' to indicate how to slice read
    """

    # find the original read seq id
    reads = [r for r in reads if 'window' not in r.seqid and
            r.seqid in extension_dict['read'].seqid]

    # we have a forward and reverse read here, so we need to make sure we use the correct one. This can be done
    # by comparing the first characters of the seqids with the extension.
    read = [r for r in reads if r.seqid[0] == extension_dict['read'].seqid[0]]

    read = read[0].read
    if side == 'left':
        read = read[:-extension_dict['read'].read_length]
    elif side == 'right':
        read = read[extension_dict['extension']:]

    print(read)