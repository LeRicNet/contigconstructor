"""functions for reading data"""

# Load query file in and convert to a simple image.
def read_query_file(query_file_path, return_all=False):
    with open(query_file_path, "r") as f:
        query = f.readlines()
    if not return_all:
        return query[-1]
    else:
        return query

def hash_fn(string, base=4, encoding={'A': 0, 'C': 1, 'G': 2, 'T': 3}):
    """simple modifiable hashing function"""
    k = 0
    for i in range(len(string)):
        k += encoding[string[i]] * base ** i
    return k

def logmessage(value, verbose):
    """simple wrapper for message handling
    value (str): standard formatted python print string
    verbose (bool): whether to print or not
    """
    if verbose:
        print(value)