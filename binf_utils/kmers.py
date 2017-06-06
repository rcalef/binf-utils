#! /usr/bin/env python2.7
"""Functions for generating kmers from sequences or an alphabet.


This module contains four functions intended for use in constructing an order k
stochastic Markov model for sequences of characters. In order to estimate 
parameters for an order k Markov model in a simple manner, one typically wants
to count the occurences of k-mers over some training set of sequences. 
The functions provided in this module can be used to gather this data, 
and are listed below:

count_kmers:         Takes in a sequence of characters, and returns a Counter
                     object containing counts of occurences of all kmers in 
                     the sequence.

kmers_from_sequence: Splits a sequence into kmers and yields each kmer in the
                     sequence, one at a time.

generate_all_kmers:  Takes in an alphabet, a start character, and a stop
                     character, and yields all valid kmers, one at a time.

check_alphabet:      Checks an alphabet, and start and stop characters to make
                     sure they are valid.

For more information on each of the functions, see their individual docstrings.
"""
from __future__ import print_function
import string
from collections import Counter
from itertools import product

def count_kmers(sequence, k, start = "^", stop = "$", ignore_case=True):
    """Counts kmers in a sequence, optionally adding start and stop characters.

    Inputs:
        sequence    - A string containing the sequence to be split into kmers,
                      without any start or stop characters.
                      
        k           - The length of kmers to be counted, i.e. k = 1 will count 
                      occurences of individual characters,

        start       - The character used to represent the start of the sequence.
                      Defaults to '^'. If set to None, then no start character
                      is used.

        stop        - A character representing the end of the sequence.
                      Defaults to '$'. If set to None, then no stop character
                      is used.

        ignore_case - A Boolean specifying whether or not to preserve case of 
                      letters in 'sequence', defaults to True. If set to False,
                      the 3-mers "AAA" and "AaA" will be treated as two distinct
                      3-mers, if True, all kmers will be converted to uppercase.
    Output:
        counts   - A Counter object where keys are kmers and values are the
                   number of occurences of that kmer.
    
    count_kmers will split 'sequence' in to kmers of length 'k', prepending and 
    appending the start and stop characters as necessary, and the occurences of
    each kmer will be counted. 
    """
    counts = Counter()
    if ignore_case:
        sequence = sequence.upper()

    # If using start and stop characters, then pad with the appropriate 
    # number of characters. 
    pad = k - 1 if k > 1 else 1
    if start:
        sequence = (start * pad) + sequence
    if stop:
        sequence = sequence + (stop * pad)

    for kmer in kmers_from_sequence(sequence,k):
        counts[kmer] += 1

    return counts

def kmers_from_sequence(sequence, k):
    """Generator to yield kmers from a sequence."""
    for start in range(len(sequence) - k + 1):
        yield sequence[start:start+k]

def generate_all_kmers(k, alphabet, start = "^", stop = "$"):
    """Generates all possible kmers over the given alphabet.
    
    Inputs:
        k        - An integer specifying the length of kmers to be generated.

        alphabet - Any iterable containing the alphabet of valid characters
                   for kmers, should not contain the start and stop characters.

        start    - The character used to represent the start of the sequence.
                   If None, then no start character is used.

        stop     - The character used to represent the end of the sequence,
                   If None, then no stop character is used.

    Outputs:
        kmer     - A string of length 'k' consisting of characters from the
                   specified alphabet, and the start and stop characters.
               All such valid kmers will be generated.
    
    generate_all_kmers takes in an alphabet, start and stop characters, and an 
    integer 'k', and yields all valid kmers over this alphabet, one kmer at a 
    time. A valid kmer comes from one of the following four cases (in the descriptions
    below, n is an integer such that 0<n<k). Examples are given for an alphabet of
    "ABC", '^' as the start character, '$' as the stop character, and k=4:
        1. A string of length k consisting only of characters from the alphabet.
            e.g. "AAAA", "AAAB", "AAAC" ...
        2. k-n characters from the alphabet prefixed by n start characters.
            e.g. "^AAB", "^^AA", "^^^A", BUT NOT "^^^^"
        3. k-n characters from the alphabet suffixed by n stop characters.
            e.g. "CCA$", "CA$$", "A$$$", BUT NOT "$$$$" 
        4. k-n characters from the alphabet prefixed by x start characters
           and suffixed by y stop characters, where x + y = n, and x and y are
           both at least 1. For this case only, n is allowed to equal k, so as 
           to generate kmers corresponding to an empty sequence.
                e.g. "^AA$", "^AB$", "^^^$", "^^$$", "^$$$"
    """
    # Generate all kmers for case 1, if not using start and stop characters,
    # then we're done.
    for kmer in product(alphabet, repeat=k):
        yield ''.join(kmer)

    # If using a start character, handle case 2
    if start:
        for num in range(1,k):
            for sub_mer in product(alphabet, repeat=k-num):
                sub_mer = ''.join(sub_mer)
                kmer = (num * start) + sub_mer
                yield kmer

    # If using a stop character, handle case 3
    if stop:
        for num in range(1,k):
            for sub_mer in product(alphabet, repeat=k-num):
                sub_mer = ''.join(sub_mer)
                kmer = sub_mer + (num * stop)
                yield kmer

    # If using both start and stop characters, handle case 4
    if start and stop:
        for num_starts in range(1,k):
            for num_stops in range(1,k - num + 1): 
                for middle_seq in product(alphabet, 
                                          repeat = k -  num_starts - num_stops):
                    kmer = (start * num_starts) + \
                           ''.join(middle_seq) + \
                           (stop * num_stops)
                    yield kmer

class AlphabetException(Exception):
    pass

def check_alphabet(alphabet,start,stop):
    """Checks alphabet for invalid characters.

    Inputs:
            alphabet - Any iterable supporting the 'in' operator and containing the
                       alphabet of allowed characters for FASTA sequences.

            start    - The character used to represent the start of a FASTA sequence.

            stop     - The character used to represent the end of a FASTA sequence.

    Output:
            none, raises exception if problems with the alphabet are found.
    """
    if start in alphabet or stop in alphabet:
        raise AlphabetException("Alphabet contains start and/or "
                                "stop characters.")

    # Allowed characters are all printable, non-whitespace characters,
    # check that the alphabet is contained in this set.
    allowed_chars = set(string.printable)-set(string.whitespace)
    if not set(alphabet) <= allowed_chars:
        raise AlphabetException("Alphabet contains either whitespace "
                                "characters, or non-printable characters.")
