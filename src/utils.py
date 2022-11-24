"""These are utility methods to help with sequence records and string operation

"""
from pathlib import Path
from typing import Dict

from Bio import SeqIO

CODON_LENGTH = 3
START_CODONS = ("ATG",)
STOP_CODONS = ("TAA", "TAG", "TGA")


def convert_index_pair_to_position_pair(index_pair):
    r"""Convert a pair of Python indices to pair positions

    Parameters
    ----------
    index_pair : tuple[int, int]
        Pair of Python indices

    Returns
    -------
    tuple[int, int]
        Parameter values

    >>> convert_index_pair_to_position_pair((0, 10))
    (1, 11)
    """
    return tuple(val + 1 for val in index_pair)


def key_by_sequence_length(seq_rec):
    r"""A key function that returns the length of the sequence record

    Parameters
    ----------
    seq_rec : SeqIO.SeqRecord

    Returns
    -------
    int
        Length of the sequence
    """
    return len(seq_rec)


def get_sequence_records(filename, fmt):
    r"""Get the sequence records from a file

    Parameters
    ----------
    filename : str
        The path to the input file
    fmt : str
        The format of the file as supported by Bio.SeqIO

    Returns
    -------
    list[SeqIO.SeqRecord]
        All sequences read from the file

    Raises
    ------
    IOError
        If the file does not exist
    """
    with Path(filename).open(encoding="utf-8") as sequence_file:
        return list(SeqIO.parse(sequence_file, fmt))


def get_forward_rf_table(seq_str, num):
    r"""Build a forward reading frame (RF) table codon index list from a DNA sequence.
    **Gaps are ignored in this treatment.**

    To use this with a sequence string to find open reading frame, use successful values of the returned list like

    seq_str[out[0] : out[1]] to get the first codon

    Parameters
    ----------
    seq_str : str
        Sequence string
    num : int

    Returns
    -------
    list[int]
        A codon reading frame mapper

    Raises
    ------
    ValueError
        When the reading frame value `num` is not in (1, 2, 3)
    """
    if num == 1:
        rf_table = list(range(0, len(seq_str) + 1, 3))
        if len(seq_str) % 3 != 0:
            rf_table.append(len(seq_str))
    elif num in (2, 3):
        rf_table = [0] + list(range(num - 1, len(seq_str) + 1, 3))
        if len(seq_str) not in rf_table:
            rf_table.append(len(seq_str))
    else:
        raise ValueError(
            "Only the forward reading frames are supported for this method"
        )
    return rf_table


def get_forward_rf_tables(rec):
    r"""Get a mapper between reader frame (RF) number and codon index list

    Parameters
    ----------
    rec : SeqIO.SeqRecord
        Sequence record

    Returns
    -------
    Dict[int, list[int]]
        A codon reading frame map with keys (1, 2, 3)
    """
    return {num: get_forward_rf_table(str(rec.seq), num) for num in range(1, 4)}


def get_all_forward_orf(rec, rf_tables):
    r"""Get all possible forward open reading frame (ORF) given a reading frame (RF) mapper

    Parameters
    ----------
    rec : SeqIO.SeqRecord
        Sequence record
    rf_tables : dict[int, list[int]]
        A codon reading frame mapper

    Returns
    -------
    Dict[int, list[tuple[int, int]]]
        Map between RF number and possible ORF
    """
    orf_by_num = dict.fromkeys(rf_tables.keys())
    for rf_num, rf_table in rf_tables.items():
        start_codon_index_pairs = []
        stop_codon_index_pairs = []
        # Get all start and stop codon positions
        for rf_table_index in range(len(rf_table) - 1):
            codon_index_start = rf_table[rf_table_index]
            codon_index_end = rf_table[rf_table_index + 1]
            test_codon = rec.seq[codon_index_start:codon_index_end]
            if test_codon in START_CODONS:
                start_codon_index_pairs.append((codon_index_start, codon_index_end))
            elif test_codon in STOP_CODONS:
                stop_codon_index_pairs.append((codon_index_start, codon_index_end))
        # print("start codon index pairs = ", start_codon_index_pairs)
        # print("stop codon index pairs = ", stop_codon_index_pairs)
        # Remove stop codons that precede start codons
        if len(start_codon_index_pairs) > 0:
            while (
                len(stop_codon_index_pairs) > 0
                and stop_codon_index_pairs[0][0] < start_codon_index_pairs[0][0]
            ):
                stop_codon_index_pairs.pop(0)

        orf = []
        # Pair-wise match for possible ORF
        while len(stop_codon_index_pairs) > 0 and len(start_codon_index_pairs) > 0:
            start_index_pair = start_codon_index_pairs[0]
            stop_index_pair = stop_codon_index_pairs[0]
            if start_index_pair[0] < stop_index_pair[0]:
                orf.append((start_index_pair[0], stop_index_pair[1]))
                start_codon_index_pairs.pop(0)
            else:
                stop_codon_index_pairs.pop(0)
        orf_by_num[rf_num] = orf
    return orf_by_num


def get_all_subsequences(rec, length):
    r"""Get the unique set of subsequences of fixed length

    Parameters
    ----------
    rec : SeqIO.SeqRecord
        Sequence record
    length : int
        Length of sequence

    Returns
    -------
    set[str]
        Set of all subsequences
    """
    return set(
        [str(rec.seq)[index : index + length] for index in range(len(rec) - length + 1)]
    )


def count_overlapping_occurrences(seq, subseq):
    r"""Find all overlapping occurrences of a subsequence. The algorithm taken from Geeks For Geeks URL
    https://www.geeksforgeeks.org/python-count-overlapping-substring-in-a-given-string/

    Parameters
    ----------
    seq : str
        Sequence string
    subseq : str
        Subsequence string

    Returns
    -------
    tuple[int, list[tuple[int, int]]]
        Occurrence count and occurrence indices
    """
    count = 0
    start = 0
    occur_indices = []
    if len(subseq) > 0 and len(seq) > 0:
        # Search through the string till
        # we reach the end of it
        while start < len(seq):

            # Check if a substring is present from
            # 'start' position till the end
            pos = seq.find(subseq, start)

            if pos != -1:
                occur_indices.append((pos, pos + len(subseq)))
                # If a substring is present, move 'start' to
                # the next position from start of the substring
                start = pos + 1

                # Increment the count
                count += 1
            else:
                # If no further substring is present
                break

    assert count == len(occur_indices)
    return count, occur_indices
