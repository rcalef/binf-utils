#! /usr/bin/env python3
"""This module contains two classes which can be used to calculate optimal local or global pairwise sequence alignments using an affine gap score.

This module contains two classes intended for calculating an optimal global or local 
pairwise sequence alignment using a provided substitution matrix and affine gap cost
parameters. The classes provided in this module are listed below:

LocalAligner  - A class defining an object that can be used for calculating optimal
                 local alignment scores, as well as outputting the alignment in A2M format.
                 Also can be used to score a local alignment in A2M format using the
                 given parameters.

GlobalAligner - A class defining an object that can be used for calculating optimal
                 global alignment scores, as well as outputting the alignment in A2M 
                 format. Also can be used to score a global alignment in A2M format using
                 the given parameters.

For a much more detailed description of the two classes, and their member attributes and
functions, please see each class' respective docstring.
"""
from __future__ import print_function
import string
import sys


from .exceptions import ParsingException,AlignmentException
from .utils import max_argmax, open_file
from .constants import neg_inf

class BaseAffineAligner:
    """Base class for calculating optimal alignments with affine gap score.

    The BaseAffineAligner class defines an interface for defining objects
    for calculating the optimal alignment of two sequences an affine gap
    score under the given set of parameters. 

    Each instance of a BaseAffineAligner subclass must be initialized with a
    substitution matrix in the format used by matblas and BLAST:
        #starts with any number of comments beginning with the '#' character.
        ...
        #more comments
            A  B  C ... X
        A   4 -2  1    -4
        B  -2  4 -1    -9 
        ...
        X  -4 -9  1     4

    The key points are:
        1. The first non-comment line is the row header containing a 
           white-space delimited list of the legal characters.
        2. Each line after the row header is a row of the substitution matrix
           itself. These lines begin with a legal character Y followed by 
           white-space delimited substitution scores for the pair where the 
           i'th score is the substitution score for (Y,C_i) where C_i is the 
           i'th character in the row header.

    Sequences passed to a BaseAffineAligner subclass should contain only legal
    characters as defined by the provided substitution matrix.

    Each instance of the BaseAffineAligner class contains nine attributes:

        self._sub_matrix - A dict storing the provided substitution matrix 
                           with tuples of characters as keys. 

        self._gap_open   - An integer specifying the open gap penalty.

        self._gap_extend - An integer specifying the gap extension penalty.

        self._double_gap - An integer specifying the double gap penalty.

        self.alphabet    - A list containing the legal characters for 
                           sequences to be aligned. Obtained from the 
                           substitution matrix row header.

        self._alignment  - A dict storing the currently cached alignment matrix.

                           Empty dict if the align() function has not been 
                           called yet. Otherwise self._alignment[(i,j)] will
                           be a tuple of three values 
                           (M_score, Ir_score, Ic_score):

                             M_score:  An integer giving the best score of an
                                       alignment of X_l...X_i and Y_k...Y_j 
                                       ending with X_i aligned to Y_j, where 
                                       X_i is the i'th character of X, and 
                                       l <= i and k <= j are the beginning
                                       indices of the alignment.

                             Ir_score: An integer giving the best score of an
                                       alignment of X_l...X_i and Y_k...Y_j 
                                       ending with X_i aligned to a gap.
                                       
                             Ic_score: An integer giving the best score of an
                                       alignment of X_l...X_i and Y_k...Y_j
                                       ending with Y_j aligned to a gap.

                            The score of the optimal alignment will be stored
                            in
                              self._alignment[(len(X),len(Y))]

                            If the align() method is called with either 
                            sequence being an empty sequence, then the 
                            self._alignment will be an empty dict.
        
        self._traceback - A dict storing the currently cached traceback matrix.

                          Empty if the align() function has not been called
                          yet. Otherwise self._traceback[(i,j)] will be a tuple
                          of three strings (M_source, Ir_source, Ic_source),
                          each of which can have one of four values: 
                          "M", "Ir", "Ic", "Start". The value specifies the
                          source of the score for the corresponding value in
                          self._alignment[(i,j)]. 

                          self._traceback[(len(X),len(Y)] is a special cell
                          that is used to store the position at which the
                          optimal alignment begins, and the source from which
                          the optimal alignment was reached:
                            self._traceback[(len(X),len(Y)] = (best_row,
                                                               best_col,
                                                               best_source)
                            e.g. (12,11,"M")

       self.row_seq     - The currently cached sequence being used as the
                          master sequence in the alignment. None if no
                          alignment has been calculated yet.
 
       self.row_seq     - The currently cached sequence being used as the
                          slave sequence in the alignment. None if no
                          alignment has been calculated yet.

    The BaseAffineAligner class also defines a constructor and two methods:

        BaseAffineAligner     - An initializer method that requires a path to
                                the file containing the desired substitution
                                matrix in the format described above.
                                Optionally takes in integers to be used for the
                                gap open, gap extension, and double gap 
                                penalties.

        read_subst_matrix     - Takes in a path to a file containing a 
                                substitution matrix in the format defined above
                                and returns (matrix,alphabet) where 'matrix' is
                                dict whose keys are tuples of two characters
                                and the value is the substitution score, and
                                'alphabet' is a set containing all the
                                characters in the substitution matrix.

        align                 - Takes in two sequences containing only 
                                characters present in the substitution matrix,
                                and calculates and stores alignment and
                                traceback matrices for these two sequences. 
                                Returns the score of the optimal alignment
                                with the given parameters.

        traceback_col_seq     - Returns the second sequence passed to the 
                                align() function in A2M format according to
                                the calculated alignment matrix. Will raise an
                                exception if called before align() has been
                                called with this instance of a
                                BaseAffineAligner subclass.

    To customize the above methods in a subclass, one can define the
    following methods that make up the interface that describe how to calculate
    the alignment:

        _init_alignment - Called at the beginning of the align() function.
                          Used to set up any boundary conditions for the
                          alignment.
        
        _get_Ic_score   
        _get_Ir_score
        _get_M_score    - These methods define how to calculate the
                          corresponding score for a given cell in the
                          alignment matrix.

        _per_cell       - Called during the align() function after each
                          time a cell is populated (excluding any cells
                          populated in _init_alignment).

        _post_alignment - Called at the end of the align() function after
                          the alignment matrix has been fully populated.
    """
    def __init__(self, subst_fn, open_penalty = 12, extend = 1, double = 3):
        """Constuct a new GlobalAligner object from a substitution matrix.

        Inputs:
            subst_fn     - Path to the file containing substitution matrix
                           data in the format described in the GlobalAligner
                           class docstring.

            open_penalty - An integer specifying the desired gap open penalty.
                           Defaults to 12.

            extend       - An integer specifying the desired gap extension
                           penalty. Defaults to 1.
            
            double       - An integer specifying the desired double gap penalty.
                           Defaults to 3.
        Output:
            _align       - A BaseAffineAligner object initialized with the
                           parameters provided as inputs.
        """
        self._gap_open   = open_penalty
        self._gap_extend = extend
        self._double_gap = double
        self._alignment  = dict()
        self._traceback  = dict()
        self._row_seq    = None
        self._col_seq    = None
        matrix, alphabet = self.read_subst_matrix(subst_fn)
        self._sub_matrix = matrix

        self.alphabet    = alphabet

    @staticmethod
    def read_subst_matrix(fn):
        """Read in a substitution matrix and return matrix and residues."""
        fh = open_file(fn)
        matrix = {}

        # Take first non-comment line as column headers
        for line in fh:
            if line.startswith("#"):
                continue
            break
        row_header = line.split()
        alphabet = set(row_header)

        # Read each row of the substitution matrix, storing score
        # based off the column headers
        for line in fh:
            row = line.split()
            first_char = row[0]
            if first_char not in alphabet:
                raise ParsingException("Substition matrix contains row that "
                                       "does not match any column: %s"
                                       % (char))

            for col_num,score in enumerate(row[1:]):
                second_char = row_header[col_num]
                matrix[(first_char,second_char)] = int(score)

        return matrix,alphabet


    def sub_score(self, charA, charB):
        """Get substitution score for characters."""
        score = self._sub_matrix.get((charA,charB))
        if score is None:
            raise AlignmentException("No substitution score for (%s,%s) in "
                                     "provided substitution matrix." 
                                     % (charA,charB))
        return score

    def reset(self):
        """Clear cached data."""
        self._alignment.clear()
        self._traceback.clear()
        self._row_seq = None
        self._col_seq = None


    def _init_alignment(self):
        raise AlignmentException("Calling align() from a BaseAffineAligner"
                                 "object!")

    def _get_Ic_score(self,row,col):
        raise AlignmentException("Calling align() from a BaseAffineAligner"
                                 "object!")

    def _get_M_score(self,row,col):
        raise AlignmentException("Calling align() from a BaseAffineAligner"
                                 "object!")

    def _get_Ir_score(self,row,col):
        raise AlignmentException("Calling align() from a BaseAffineAligner"
                                 "object!")

    def _per_cell(self,row,col):
        raise AlignmentException("Calling align() from a BaseAffineAligner"
                                 "object!")

    def _post_alignment(self,row,col):
        raise AlignmentException("Calling align() from a BaseAffineAligner"
                                 "object!")

    def align(self,row_seq,col_seq):
        """Calculate the optimal global alignment of two sequences.

        Inputs:
            row_seq    - A string containing the master sequence to align 
                         'col_seq' against.
                      
            col_seq    - A string containing the slave sequence to align to 
                         'row_seq'.
       
            WARNING: Both 'row_seq' and 'col_seq' must contain only characters
                     present in the substitution matrix used to initialize the 
                     calling GlobalAligner object or an AlignmentException will
                     be raised.

        Output:
            best_score - An integer giving the score of an optimal alignment of
                         'row_seq' and 'col_seq' under the parameters provided
                         to the calling object. 

        The align() function uses the parameters used to initialize the calling
        object to calculate the alignment matrices for an alignment of 'row_seq'
        and 'col_seq' using an affine gap penalty. The alignment matrices will
        be cached in the 'self._alignment' attribute of the calling object, 
        that is, after calling aligner.align(s1,s2), 'aligner._alignment' will
        contain the alignment matrices for the sequences s1 and s2. 

        While the alignment matrix is being calculated, a traceback matrix will
        be generated simultaneously, and stored in 'self._traceback' for later 
        use.

        This function only calculates the alignment matrices and returns the
        score of an optimal alignment. To obtain the actual sequence of the
        alignment in A2M format, use the self.traceback_col_seq() function.

        For a more rigorous description of the meaning of the elements of 
        'self._alignment', see the BaseAffineAligner class docstring.
        """

        # Clear any data from previous alignments.
        self.reset()

        self._row_seq = row_seq.upper()
        self._col_seq = col_seq.upper()
        row_len = len(self._row_seq)
        col_len = len(self._col_seq)

        if row_len == 0 or col_len == 0:
            self._traceback[(row_len,col_len)] = (-1,-1)
            return 0

        self._init_alignment()


        # Populate the alignment matrix
        for row in range(row_len):
            for col in range(col_len):
                # For each cell, calculate the three alignment scores 
                # ("M","Ir","Ic") corresponding to an alignment ending with
                # X_i,Y_j aligned, X_i aligned to a gap, or Y_j aligned to a gap
                # respectively.

                Ic_score, Ic_source = self._get_Ic_score(row,col)
                M_score,M_source = self._get_M_score(row,col)
                Ir_score, Ir_source = self._get_Ir_score(row,col)

                self._alignment[(row,col)] = (M_score,Ir_score,Ic_score)
                self._traceback[(row,col)] = (M_source,Ir_source,Ic_source)

                self._per_cell(row,col)

        # Finally, we store the index of the cell containing the best score,
        # or (-1,-1) if no local alignment with a score greater than zero.
        self._post_alignment()
        return self._alignment[(row_len,col_len)]

    def traceback_col_seq(self):
        """Get the alignment of the cached slave (column) sequence in A2M format

        Inputs:
             None

        Outputs:
             col_a2m - A string giving the slave sequence (the second sequence
                       passed to the align() function) in A2M format. Describes
                       an optimal local alignment of the slave sequence to the
                       master sequence (the first sequence passed to the align()
                       function).

             WARNING: This function should only be called after the align()
                      function has been called from the same 
                      object, otherwise there's no cached alignment to 
                      traceback on. If this function is called before the 
                      align() function has been called, an AlignmentException
                      will be raised.

        traceback_col_seq() is used to perform the traceback step of the
        alignment algorithm, returning the slave sequence in A2M format,
        providing a description of an optimal local aligment of the slave
        sequence to the master sequence.
        """
        # If no cached alignment, throw exception
        if self._row_seq is None or self._col_seq is None:
            raise AlignmentException("traceback_col_seq called before the "
                                     "align() function.")

        row_len = len(self._row_seq)
        col_len = len(self._col_seq)
        col_a2m = []

        # curr_row,curr_col will be -1 if best alignment is the
        # empty alignment, first source must be a match if not
        curr_row, curr_col, curr_source = self._traceback[(row_len,col_len)]

        # Building up the aligned sequence from the end, so add any trailing 
        # gaps or insertions relative to the master sequence
        col_a2m.extend(['-' for i in range(curr_row + 1, row_len)])
        col_a2m.extend(self._col_seq[curr_col + 1:col_len][::-1].lower())

        while(True):
            # If either index is less than 0, than we've run off the end of
            # the alignment, so break out of the traceback and add any
            # remaining characters.
            if curr_row < 0 or curr_col < 0:
                break

            M_source, Ir_source, Ic_source = self._traceback[(curr_row,curr_col)]
            seq_char = self._col_seq[curr_col]
            # If source is "Start", then this is the end of the alignment
            if curr_source == "Start":
                break

            # Match, so add character to traceback and decrement row
            # and column
            elif curr_source == "M":
                next_char = seq_char
                next_source = M_source
                curr_row -=1
                curr_col -=1

            # Gap relative to master
            elif curr_source == "Ir":
                next_char = "-"
                next_source = Ir_source
                curr_row -= 1

            # Insertion relative to master
            elif curr_source == "Ic":
                next_char = seq_char.lower()
                next_source = Ic_source
                curr_col -= 1

            col_a2m.append(next_char)
            curr_source = next_source

        # Add in any gaps or characters from the slave sequence and return
        # reverse since we're building up the alignment from the end
        col_a2m.extend(self._col_seq[0:curr_col + 1][::-1].lower())
        col_a2m.extend(['-' for i in range(curr_row + 1)])
        return ''.join(col_a2m[::-1])

class GlobalAligner(BaseAffineAligner):
    """Calculates optimal global alignments with an affine gap score.

    The GlobalAligner class defines an object used for calculating the optimal
    global alignment of two sequences with the given set of parameters. 
    """
    def __init__(self, subst_fn, open_penalty = 12, extend = 1, double = 3):
        """Constuct a new GlobalAligner object from a substitution matrix.
        """
        super().__init__(subst_fn, open_penalty, extend, double)


    def _init_alignment(self):
        """Set boundary conditions for global alignment."""
        # For global alignment, we can align the beginning of either sequence
        # to gaps, so we use a "-1" row and column to indicate alignment to a
        # gap. The (-1,-1) cell can't extend a previous gap, so we set Ir and
        # Ic values to neg_inf, and M to 0. For the (i,-1) and (-1,j) cells,
        # only the Ir and Ic scores are valid respectively, corresponding to
        # opening and extending a gap in either the row or column sequence.
        #
        # No need to set traceback values, as we'll exit traceback when we
        # reach any of these "-1" cells.
        self._alignment[(-1,-1)] = (0, neg_inf, neg_inf)
        for i in range(len(self._row_seq)):
            self._alignment[(i,-1)] = (neg_inf,
                                       -(self._gap_open + i * self._gap_extend),
                                       neg_inf)
        for j in range(len(self._col_seq)):
            self._alignment[(-1,j)] = (neg_inf,
                                       neg_inf,
                                       -(self._gap_open + j * self._gap_extend))

    def _get_Ic_score(self,row,col):
        """Get score for Ir matrix in a given alignment cell."""
        # For the Ic score, we have three possiblities for the new
        # score:
        # opening a gap, extending a gap, or switching to a gap in the
        # other sequence.
        M_up, Ir_up, Ic_up = self._alignment[(row,col-1)]
        return max_argmax(M  = M_up - self._gap_open,
                          Ic = Ic_up - self._gap_extend,
                          Ir = Ir_up - self._double_gap)

    def _get_M_score(self,row,col):
        """Get score for M matrix in a given alignment cell."""
        # For the M score, all possible scores correspond to aligning
        # X_i and Y_j, so just take a max of the possible ways to 
        # reach this alignment:
        #   previous  match 
        #   previous gap in column sequence
        #   previous gap in row sequence
        M_diag, Ir_diag, Ic_diag = self._alignment[(row-1,col-1)]
        M_score, M_source = max_argmax(M  = M_diag,
                                       Ic = Ic_diag,
                                       Ir = Ir_diag)
        M_score += self.sub_score(self._row_seq[row],
                                  self._col_seq[col])
        return M_score,M_source

    def _get_Ir_score(self,row,col):
        """Get score for Ir matrix in a given alignment cell."""
        # The Ir score is calculated in the same manner as the Ic 
        # score.
        M_left,Ir_left,Ic_left = self._alignment[(row-1,col)]
        return max_argmax(M  = M_left - self._gap_open,
                          Ic = Ic_left - self._double_gap,
                          Ir = Ir_left - self._gap_extend)

    def _per_cell(self,row,col):
        pass

    def _post_alignment(self):
        # Since we alignments can end with X_i and/or Y_j aligned to a 
        # gap, we have to take the max of the three final global alignment
        # scores and store the score and source. By definition, the optimal
        # global alignment ends at (i,j).
        row_len = len(self._row_seq)
        col_len = len(self._col_seq)
        final_M, final_Ir, final_Ic = self._alignment[(row_len - 1,
                                                       col_len - 1)]
        final_score, final_source = max_argmax(M  = final_M,
                                               Ic = final_Ic,
                                               Ir = final_Ir)

        self._alignment[(row_len,col_len)] = final_score
        self._traceback[(row_len,col_len)] = (row_len - 1,
                                              col_len - 1,
                                              final_source)

class LocalAligner(BaseAffineAligner):
    """Calculates optimal local alignments with an affine gap score.

    The LocalAligner class defines an object used for calculating the optimal
    local alignment of two sequences with the given set of parameters. 

    Additions to the BaseAffineAligner class are as follows:

      self._best_index - Used to store the index of the highest scoring local
                         alignment.

      self._best_score - Used to store the score of an optimal local 
                         alignment.
    """
    def __init__(self, subst_fn, open_penalty = 12, extend = 1, double = 3):
        """Constuct a new LocalAligner object from a substitution matrix."""
        super().__init__(subst_fn, open_penalty, extend, double)
        self._best_index = (-1,-1)
        self._best_score = 0

    def _init_alignment(self):
        self._best_index = (-1,-1)
        self._best_score = 0

    def _get_Ic_score(self,row,col):
        # Local alignments can't start with gaps by definition, so all
        # starting cells have Ic values of negative infinity
        if row == 0 or col == 0:
            return neg_inf,None

        # For the Ic score, we have three possiblities for the new
        # score:
        # opening a gap, extending a gap, or switching to a gap in the
        # other sequence.
        M_up, Ir_up, Ic_up = self._alignment[(row,col-1)]
        return max_argmax(M  = M_up - self._gap_open,
                          Ic = Ic_up - self._gap_extend,
                          Ir = Ir_up - self._double_gap)

    def _get_M_score(self,row,col):
        # If beginning cell, then this must be the beginning of an alignment
        # so just return the substitution score.
        if row == 0 or col == 0:
            return self.sub_score(self._row_seq[row],self._col_seq[col]),None

        # For the M score, all possible scores correspond to aligning
        # X_i and Y_j, so just take a max of the possible ways to 
        # reach this alignment:
        #   previous  match 
        #   previous gap in column sequence
        #   previous gap in row sequence
        #   starting a new alignment
        M_diag, Ir_diag, Ic_diag = self._alignment[(row-1,col-1)]
        M_score, M_source = max_argmax(M  = M_diag,
                                       Ic = Ic_diag,
                                       Ir = Ir_diag,
                                       Start = 0)
        M_score += self.sub_score(self._row_seq[row],
                                  self._col_seq[col])
        return M_score,M_source

    def _get_Ir_score(self,row,col):
        # The Ir score is calculated in the same manner as the Ic 
        # score.
        if row == 0 or col == 0:
            return neg_inf,None

        M_left,Ir_left,Ic_left = self._alignment[(row-1,col)]
        return max_argmax(M  = M_left - self._gap_open,
                          Ic = Ic_left - self._double_gap,
                          Ir = Ir_left - self._gap_extend)

    def _per_cell(self,row,col):
        # Want to keep track of index and score of best alignment seen
        # so far for easy traceback. Local alignment must end in a match
        # by definition, so only concerned with the M score.
        M_score, _, _ = self._alignment[(row,col)]
        if M_score > self._best_score:
            self._best_score = M_score
            self._best_index = (row,col)

    def _post_alignment(self):
        # Need to store index and score of best local alignment
        row_len = len(self._row_seq)
        col_len = len(self._col_seq)
        
        self._alignment[(row_len,col_len)] = self._best_score
        self._traceback[(row_len,col_len)] = self._best_index + ("M",)
