#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This module is designed for analyzing RNA secondary structures. It includes functionality to
identify base pairing positions in RNA structures and calculate statistical scores like MCC (Matthews 
correlation coefficient), precision, recall, and F-score for comparing two RNA secondary structures.

The module defines a class `SimScore` for calculating these statistical measures, and a utility
function `pairing_positions` to map the base pairs in RNA secondary structures.
"""

import math

OPEN_BRACKETS = {"(": "0", "[": "1", "<": "2", "{": "3", "A": "4", "B": "5", "C": "6", "D": "7", "E": "8"}
CLOSE_BRACKETS = {")": "0", "]": "1", ">": "2", "}": "3", "a": "4", "b": "5", "c": "6", "d": "7", "e": "8"}


def pairing_positions(s1):
    """
    Maps the base pairing positions in an RNA secondary structure.

    Parameters:
    s1 (str): An RNA secondary structure represented as a string with brackets for base pairs.

    Returns:
    dict: A dictionary with key-value pairs representing the base pairing positions.
    """
    final_pairs = []
    l_opens = []
    l_closes = []
    for i in range(len(s1)):
        if s1[i] in OPEN_BRACKETS:
            l_opens.append([i, OPEN_BRACKETS[s1[i]]])
        if s1[i] in CLOSE_BRACKETS:
            l_closes.append([i, CLOSE_BRACKETS[s1[i]]])
        if s1[i] == '.' or s1[i] == '-':
            final_pairs.append((i, -1))

    lengh = (len(l_opens) - 1)
    for i in range(lengh, -1, -1):
        for j in range(lengh + 1):
            if (l_opens[i][-1]) == (l_closes[j][-1]) and (l_closes[j][0] > l_opens[i][0]):
                final_pairs.append((l_opens[i][0], l_closes[j][0]))
                final_pairs.append((l_closes[j][0], l_opens[i][0]))
                del l_opens[i][-1]
                del l_closes[j][-1]

    return (dict(sorted(final_pairs)))


class SimScore:
    """
    A class for scoring similarity between two RNA secondary structures. It calculates various
    statistical scores such as MCC, precision, recall, and F-score.

    Attributes:
    ref_ss (str): The reference RNA secondary structure.
    query_ss (str): The query RNA secondary structure to compare against the reference.
    """

    def __init__(self, ref_ss, query_ss):
        """
        Initializes the SimScore instance with reference and query RNA secondary structures.

        Parameters:
        ref_ss (str): The reference RNA secondary structure.
        query_ss (str): The query RNA secondary structure.
        """
        self.ref_ss = ref_ss
        self.query_ss = query_ss

    def find_basepairs(self):
        """
        Identifies and stores the base pairing positions in both the reference and query RNA secondary structures.
        """
        self.bp_dict_r = pairing_positions(self.ref_ss)
        self.bp_dict_q = pairing_positions(self.query_ss)

    def cofusion_matrix(self):
        """
        Calculates the confusion matrix for the comparison between the reference and query RNA secondary structures.

        The confusion matrix includes true positives (tp), false positives (fp), true negatives (tn), and false negatives (fn).
        """
        tp = 0
        fp = 0
        tn = 0
        fn = 0
        for i in range(len(self.bp_dict_r)):
            if self.bp_dict_r[i] == self.bp_dict_q[i] and self.bp_dict_r[i] != -1:
                tp += 1
            if self.bp_dict_r[i] == self.bp_dict_q[i] and self.bp_dict_r[i] == -1:
                tn += 1
            if self.bp_dict_r[i] != self.bp_dict_q[i]:
                if self.bp_dict_r[i] == -1:
                    fp += 1
                else:
                    fn += 1
        self.conf_mat = (tp, fp, fn, tn)

    def mcc(self):
        """
        Calculates the Matthews correlation coefficient (MCC) based on the confusion matrix.

        Returns:
        float: The MCC value, rounded to three decimal places.
        """
        tp, fp, fn, tn = self.conf_mat
        if (tp == 0 and fp == 0 and fn == 0 and tn != 0):
            numerator = 1
            denominator = 1
        else:
            numerator = ((tp * tn) - (fp * fn))
            denominator = math.sqrt((tp + fp) * (tp + fn) * (tn + fn) * (tn + fp))
        epsilon = 0.00001
        return round(numerator / (denominator + epsilon), 3)

    def recall(self):
        """
        Calculates the recall score based on the confusion matrix.

        Returns:
        float: The recall value, rounded to three decimal places.
        """
        tp, fp, fn, tn = self.conf_mat
        return round((tp / (tp + fn + 0.001)), 3)

    def precision(self):
        """
        Calculates the precision score based on the confusion matrix.

        Returns:
        float: The precision value, rounded to three decimal places.
        """
        tp, fp, fn, tn = self.conf_mat
        return round((tp / (tp + fp + 0.001)), 3)

    def fscore(self):
        """
        Calculates the F-score based on the precision and recall scores.

        Returns:
        float: The F-score, rounded to four decimal places.
        """
        makhraj = self.precision() + self.recall()
        if makhraj < 0.001:
            makhraj = 0.001

        return round(2 * (self.precision() * self.recall() / (makhraj)), 4)

    def mcc_reverse(self):
        """
        Calculates the negative value of the Matthews correlation coefficient (MCC).

        Returns:
        float: The negative MCC value.
        """
        return -self.mcc()

    def recall_reverse(self):
        """
        Calculates the negative value of the recall score.

        Returns:
        float: The negative recall value.
        """
        return -self.recall()

    def precision_reverse(self):
        """
        Calculates the negative value of the precision score.

        Returns:
        float: The negative precision value.
        """
        return -self.precision()

    def fscore_reverse(self):
        """
        Calculates the negative value of the F-score.

        Returns:
        float: The negative F-score.
        """
        return -self.fscore()
