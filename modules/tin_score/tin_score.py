###############################################################################
#
#   Main script which does the job of parsing values and plotting graphs
#   for tin scores
#
#   AUTHOR: Krish Agarwal
#   AFFILIATION: University_of_Basel
#   CONTACT: akrish136@gmail.com
#   CREATED: 1-11-2020
#   LICENSE: Apache_2.0
#
###############################################################################

from __future__ import print_function
import logging
import os

from multiqc.plots import linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    This class is instantiated in setup.py file and contains
    all the functions required by the tin-score plugin.
    """

    def __init__(self):
        """
        Initializes the parent function of MultiQC and calls
        relevant functions for the searching of files and the
        plotting of graphs.
        """
        super(MultiqcModule, self).__init__(
            name="tin-score",
            anchor="tin-score",
            href="",  # to add
            info="",  # to add
        )

        self.samples = []
        self.number = 0
        self.findLogs()

    def get_sample_name(self, f: dict):
        """
        Extracts the directory from which the file is being parsed.
        """
        sample = os.path.split(os.path.split(f["root"])[0])[1]
        return sample

    def findLogs(self):
        """
        Find files matching with the regex in and passes them to the
        parsing function one by one and also plots the graph in the end.
        """
        data = dict()

        for f in self.find_log_files("tin-score"):
            sample = self.get_sample_name(f)
            data.setdefault(sample, {})
            if sample not in self.samples:
                self.number = self.number + 1
                self.samples = self.samples + [sample]
                sample_data = self.parse_tin_score_logs(f)

                data[sample].update(sample_data)

        if self.number == 0:
            raise UserWarning

        self.add_section(
            name="TIN-score",
            anchor="tin-score",
            plot=linegraph.plot(data),
        )

    def parse_tin_score_logs(self, f: dict):
        """
        Parses the tsv file by storing the count of TIN-score in
        a dictionary.
        """
        sample_data = dict()
        word = []
        words = []
        listToStr1 = []
        listToStr = []

        for char in f["f"]:
            if char != "\t" and char != "\n":
                word.append(char)
            else:
                words.append(word)
                word = []

        # Concatinating each character
        for k in words:
            listToStr1.append("".join([str(elem) for elem in k]))
        # Conactinating each word
        for k in listToStr1:
            listToStr.append("".join([str(elem) for elem in k]))

        for i in range(len(listToStr)):
            if i < 2:
                continue
            if i % 2 == 1:
                sample_data.setdefault(round(float(listToStr[i])), 0)
                sample_data[round(float(listToStr[i]))] = (
                    sample_data[round(float(listToStr[i]))] + 1
                )

        return sample_data