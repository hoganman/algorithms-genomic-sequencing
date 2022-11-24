"""The SequenceRecordList class is a collection of SeqRecords from BioPython"""
import bisect
from collections import UserList
from typing import List

from Bio import SeqIO

from . import utils


class SequenceRecordList(UserList):
    r"""Tasks for a list of sequence records. Each entry is a SeqIO.SeqRecord

    The input list is sorted by sequence length from shorted [0] to longest [-1].
    """

    def __init__(self, data=None):
        r"""Initialize the sequence record list

        Parameters
        ----------
        data : list[SeqIO.SeqRecord] | None
            The list of sequence records
        """
        self.data: List[SeqIO.SeqRecord]
        super().__init__(data)
        self.data.sort(key=lambda record: len(record.seq))
        self.sorted_ids = {
            self.data[index].id: index for index in range(len(self.data))
        }

    @classmethod
    def read_file(cls, filename, fmt):
        r"""Get the sequence records from a file

        Parameters
        ----------
        filename : str
            The path to the input file
        fmt : str
            The format of the file as supported by Bio.SeqIO
        """
        sequences = utils.get_sequence_records(filename, fmt)
        return cls(sequences)

    def get_records_of_length(self, length):
        r"""Find all the sequence records that are of the given length

        Parameters
        ----------
        length : int
            The target length of the sequence

        Returns
        -------
        tuple[SeqIO.SeqRecord]
            The found sequences that are of target length
        """
        # For really long sequences, this saves a lot of time
        start_index = bisect.bisect_left(
            self.data, length, key=utils.key_by_sequence_length
        )
        end_index = bisect.bisect_right(
            self.data, length, key=utils.key_by_sequence_length
        )
        return tuple(self.data[start_index:end_index])

    def get_longest_sequences(self):
        r"""Get the sequences of the longest length

        Returns
        -------
        tuple[SeqIO.SeqRecord]
            The found sequences that are of the longest length
        """
        target_length = len(self.data[-1])
        return self.get_records_of_length(target_length)

    def get_shortest_sequences(self):
        r"""Get the sequences of the shortest length

        Returns
        -------
        tuple[SeqIO.SeqRecord]
            The found sequences that are of the shortest length
        """
        target_length = len(self.data[0])
        return self.get_records_of_length(target_length)

    def get_by_id(self, ident):
        r"""Get the sequence record by ID

        Parameters
        ----------
        ident : str
            Sequence record ID

        Returns
        -------
        SeqIO.SeqRecord
            The found record
        """
        if ident not in self.sorted_ids:
            return None
        return self.data[self.sorted_ids[ident]]

    def get_all_forward_repeating_overlapping_subsequences_of_length(self, length):
        r"""Highly specialized method to find forward-direction, repeating, and overlapping subsequences of a fixed
        length

        Parameters
        ----------
        length : int
            The length of subsequences to find

        Returns
        -------
         dict[str, dict[str]]
            A dictionary mapping from found subsequence to a JSON-dictionary where the key and value is
            "ids" is a list of SeqRecord id, "index_pairs" is a list of tuple-pairs of the start-stop indices for the
            subsequence for the same arg of "ids", and "count" is the count of the occurrences in the entire collection.
        """
        subsequences_strings: set[str] = set()
        # Find all possible subsequences
        for seq_rec1 in self.data:
            subsequences = utils.get_all_subsequences(seq_rec1, length)
            subsequences_strings = subsequences_strings.union(subsequences)

        subsequence_counts = dict.fromkeys(subsequences_strings, None)
        # Search all records for repeating of the sequences
        for subseq in subsequence_counts.keys():
            counts_dict = {"ids": list(), "index_pairs": list(), "count": 0}
            for seq_rec2 in self.data:
                # print("record: ", seq_rec2.seq)
                # print("subseq: ", subseq)
                subseq_count, subseq_index_pairs = utils.count_overlapping_occurrences(
                    str(seq_rec2.seq), subseq
                )
                # print("pairs: ", subseq_index_pairs)
                # print("======================")
                if subseq_count == 0:
                    continue
                for index_pair in subseq_index_pairs:
                    if (
                        seq_rec2.id in counts_dict["ids"]
                        and index_pair in counts_dict["index_pairs"]
                    ):
                        continue
                    counts_dict["count"] += 1
                    counts_dict["ids"].append(seq_rec2.id)
                    counts_dict["index_pairs"].append(index_pair)
            subsequence_counts[subseq] = counts_dict
            # print(subsequence_counts)
        return subsequence_counts
