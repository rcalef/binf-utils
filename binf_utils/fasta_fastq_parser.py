#! /usr/bin/env python3.3
"""Classes for representing FASTA and FASTQ sequences and parsers.

This module defines classes for representing FASTA and FASTQ sequences, 
'FastaEntry' and 'FastqEntry' respectively. Also implemented are parsers for
FASTA files, FASTA file with an associated .qual file, and FASTQ files:
    'FastaParser', 'FastaWithQualityParser', and "FastqParser"
"""
import sys
import re
import operator
import argparse

from .utils import open_file

class ParsingException(Exception):
    pass

class FastaSeq:
    """Represents a single sequence in FASTA format.


        self.identifier - Sequence ID as a string. Sequence ID defined as 
                          everything on the header line after the '>' and 
                          up to the first whitespace character or comma. 

        self.comment    - Comment as a string, everything from the header that
                          isn't part of the ID.

        self.sequence   - Original sequence as a string. 
    """
    def __init__(self, identifier, comment, sequence):
        self.identifier = identifier
        self.comment = comment
        self.sequence = sequence

    def write_fasta(self, output, line_len = 100):
        curr_char = 0
        print(">%s %s" % (self.identifier, self.comment), file = output)
        while(curr_char < len(self.sequence)):
            print(self.sequence[curr_char:curr_char + line_len], file = output)
            curr_char += line_len

class FastqSeq(FastaSeq):
    """Represents a single FASTQ sequence.

       Contains all the same fields as 'fasta_seq' with the addition of: 

        self.scores    -  A list of integers representing the Phred quality scores 
                          for each symbol in the sequence. Specifically, there is 
                          guaranteed to be one score for each character in 
                          self.sequence, where self.scores[i] is the quality
                          score for self.sequence[i].
    """
    def __init__(self, identifier, comment, sequence, qual_scores):
        super().__init__(identifier,comment,sequence)
        self.scores=qual_scores

    def write_fastq(self, output, phred = 33):
        print("@%s %s\n%s\n+\n%s" % (self.identifier,
                                     self.comment,
                                     self.sequence,
                                     "".join(map(lambda x: chr(x+phred),
                                                 self.scores))),
                                    file = output)

    def write_qual(self, output):
        print(">%s %s\n%s" % (self.identifier, 
                              self.comment,
                              " ".join(map(str,self.scores))),
                            file = output)

class BaseFastxParser:

    @staticmethod
    def header_split(header):
        """Split a FASTA/FASTQ header line in to the ID and comment."""
        # Remove initial character and find the first whitespace or comma
        header = header.strip()[1:]
        match = re.search("[\s,]",header)

        #If no match occurs, return the header without the initial character as the sequence ID
        if match is None: 
            return(header,"")

        #Otherwise split the header about the beginning of the first match to get ID and comment
        seq_id = header[:match.start(0)]
        comment = header[match.start(0)+1:]
        return (seq_id,comment)



class FastaParser(BaseFastxParser):
    """Parser for FASTA files.

    Implemented as an iterator over the entries in the FASTA file.
    """
    fasta_allowed_bases = set(['A','B','C','D','E','F','G','H','I','K','L',
                               'M','N','P','Q','R','S','T','U','V','W','X',
                               'Y','Z','-','*'])

    def __init__(self, filename, replace_invalid = False):
        self._fasta_fh = open_file(filename)
        self.replace_invald = replace_invalid

    def __iter__(self):
        return self._get_fasta_seqs()

    def _get_fasta_seqs(self):
        line  = ""
        lines = []

        # skip to next header and get id and comment
        while not line.startswith('>'):
            line = next(self._fasta_fh)
        identifier,comment = self.header_split(line)

        for line in self._fasta_fh:
            if line.startswith('>'):
                yield FastaSeq(identifier,comment,"".join(lines))
                identifier,comment = self.header_split(line)
                lines = []
                continue

            line = line.strip()
            for pos,char in enumerate(line):
                if char.upper() not in self.fasta_allowed_bases:
                    print("WARNING Encountered symbol not "
                          "in alphabet: %s\nSequence: %s, position: %d"
                          % (char,identifier,pos),file=sys.stderr)
                    if replace_invalid:
                        line[pos] = 'N'
            lines.append(line) 

        #If we reach the end of file, then return last entry
        yield FastaSeq(identifier,comment,"".join(lines))


class FastaWithQualityParser(FastaParser):
    """Parser for FASTA file with a separate .qual file.

    Implemented as an iterator over the FASTA and .qual files, will error out
    if there's any mismatch in the number of sequences in the two files, the
    lengths of the sequences, or their identifiers.
    """
    def __init__(self, fasta_fn, qual_fn, replace_invalid = False):
        super().__init__(fasta_fn,replace_invalid)
        self._qual_fh = open_file(qual_fn)

    def __iter__(self):
        return self

    def __next__(self):
        found_fasta = True
        try:
            fasta = next(self._get_fasta_seqs())
        except StopIteration:
            found_fasta = False

        found_qual = True
        try:
            qual = next(self._get_quals())
        except StopIteration:
            found_qual = False

        if any((found_fasta,found_qual)) and not \
           all((found_fasta,found_qual)):
               raise ParsingException("ERROR: Number of FASTA sequences not "
                                      "equal to number of quality sequences.")

        if not found_fasta and not found_qual:
            raise StopIteration

        self._validate_seq_and_qual(fasta,qual)

        return FastqSeq(fasta.identifier,
                        fasta.comment,
                        fasta.sequence,
                        qual.scores)

    @staticmethod
    def _validate_seq_and_qual(fasta,qual):
        if fasta.identifier != qual.identifier:
            raise ParsingException("ERROR: FASTA sequence and associated "
                                   "quality score entry have different "
                                   "identifiers. FASTA: %s, scores: %s" 
                                   % (fasta.identifer,qual.identifier))
        if len(fasta.sequence) != len(qual.scores):
            raise ParsingException("ERROR: FASTA sequence '%s' and associated "
                                   "quality scores have different lengths. "
                                   "FASTA: %d, scores: %d" 
                                   % (fasta.identifier,
                                      len(fasta.sequence),
                                      len(qual.scores)))

    def _get_quals(self):
        in_score = False
        scores=[]
        line = ""

        # skip to next header and get id and comment
        while not line.startswith('>'):
            line = next(self._qual_fh)
        identifier,comment = self.header_split(line)

        for line in self._qual_fh:
            if line.startswith('>'):
                yield FastqSeq(identifier,comment,"",scores)

            # Read the scores, raise exception if invalid score 
            parts = line.split()
            for part in parts:
                if not part.isdigit() or int(part) < 0:
                    raise ParsingException("ERROR: Encountered invalid "
                                           "quality score for sequence %s, "
                                           "score: %s" % (identifier,part))
            scores.extend(map(int,parts))
        yield FastqSeq(identifier,comment,"",scores)

class FastqParser(BaseFastxParser):
    """ Parser for FASTQ files.

    Implemented as an iterator over the entries in the FASTQ file. Assumes
    the standard four line structure of a fastq entry: 
        1. header line beginning with '@'
        2. sequence symbols on a single line
        3. separator line beginning with '+'
        4. Phred-encoded quality scores on a single line
    """
    fastq_allowed_bases = set(['A','C','G','T','N','U','K','S','Y','M',
                               'W','R','B','D','H','V','-'])

    def __init__(self, fastq_fn, phred = 33, replace_invalid = False):
        self._fastq_fh = open_file(fastq_fn)
        self._phred = phred
        self._replace_invalid = replace_invalid

    def __iter__(self):
        return self._get_fastq_entries()

    def _get_fastq_entries(self):
        seq_id = ""
        comment = "" 
        sequence = ""
        scores = []

        while True:
            # Skip to next header line, rely on this loop raising
            # StopException to end iteration at end of file
            try:
                line = next(self._fastq_fh)
                if not line.startswith('@'):
                    raise ParsingException("ERROR: Malformed FASTQ file near "
                                           "sequence: %s" % (seq_id))
                seq_id,comment = self.header_split(line)
            except StopIteration:
                break

            try:
                sequence = next(self._fastq_fh).strip()
                for pos,symbol in enumerate(sequence):
                    if symbol not in self.fastq_allowed_bases:
                        print("WARNING Encountered symbol not "
                              "in alphabet: %s\nSequence: %s, position: %d"
                        % (char,identifier,pos),file=sys.stderr)

                        if replace_invalid:
                            line[pos] = 'N'
                _ = next(self._fastq_fh)
                score_line = next(self._fastq_fh).strip()
            except StopIteration:
                raise ParsingException("ERROR: Unexpected end of file when "
                                       "parsing FASTQ sequence: %s" % (seq_id))
            for char in score_line:
                score = ord(char) - self._phred
                if score < 0:
                    raise ParsingException("ERROR: Phred score less than zero "
                                           "(wrong Phred offset?), sequence: %s, "
                                           "character: %s" % (seq_id,char))
                scores.append(score)
            yield FastqSeq(seq_id, comment, sequence, scores)
