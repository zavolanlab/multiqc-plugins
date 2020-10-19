###############################################################################
#
#   Main script which does the job of parsing values and plotting graphs
#
#   AUTHOR: Krish Agarwal
#   AFFILIATION: University_of_Basel
#   CONTACT: akrish136@gmail.com
#   CREATED: 18-10-2020
#   LICENSE: Apache_2.0
#
###############################################################################

from __future__ import print_function
import logging
import os

from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="ALFA",
            anchor="ALFA",
            href="https://github.com/biocompibens/ALFA",
            info="is an example analysis module used for writing \
                documentation.",
        )

        self.reset()
        self.findLogs("Unique")
        self.print_alfa_charts("Unique")
        self.reset()
        self.findLogs("UniqueMultiple")
        self.print_alfa_charts("UniqueMultiple")

    def reset(self):
        """
        Resets the required variables in order to make them ready
        for the parsing other folders.
        """
        self.alfa_data = dict()
        self.alfa_data["temp"] = {}

        self.categories = dict()
        self.categories["temp"] = {}

        self.biotypes = dict()
        self.biotypes["temp"] = {}

        self.categories_temp = dict()
        self.categories_temp["temp"] = {}

        self.biotypes_temp = dict()
        self.biotypes_temp["temp"] = {}

        # Number of ALFA reports
        self.number = 0

        # To keep track of all the parsed files
        self.filesDone = []

    def findLogs(self, folder: str):
        """
        Find files matching with the regex in ALFA/custom_code.py and
        calls iterates through all these files and calls parse_alfa_logs()
        and updateDictValues() for each file.
        """
        for f in self.find_log_files("ALFA"):
            if os.path.basename(f["root"]) == folder:
                filename = self.get_filename(f)

                if filename not in self.filesDone:
                    self.number = self.number + 1
                    self.filesDone = self.filesDone + [filename]
                    self.categories_temp, self.biotypes_temp = self.parse_alfa_logs(
                        f, filename
                    )
                    self.updateDictValues(self.categories_temp, self.biotypes_temp)

                if self.number == 0:
                    raise UserWarning

    def get_filename(self, f: dict):
        """
        Extracts the main name from the filename by removing
        common parts in each filename.
        """
        fullFilename = f["fn"]
        part2 = fullFilename.find(".ALFA_feature_counts.tsv")
        part1 = fullFilename[0:part2]
        return part1

    def parse_alfa_logs(self, f: dict, filename: str):
        """
        Iterates through a file and stores biotypes' and categories'
        aggregated values in appropriate dictionaries.
        """
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

        catdict = dict()
        catdict["temp"] = {}

        biodict = dict()
        biodict["temp"] = {}

        # Adding tsv values to a dictionary
        for i_num in range(len(listToStr)):
            if i_num < 3:
                continue
            elif i_num % 3 == 0:
                cat, bio = listToStr[i_num].split(",")

                catdict.setdefault(filename, {})
                biodict.setdefault(filename, {})

                catdict[filename].setdefault(cat, 0)
                biodict[filename].setdefault(bio, 0)

                catdict[filename][cat] = catdict[filename][cat] + float(
                    listToStr[i_num + 1]
                )
                biodict[filename][bio] = biodict[filename][bio] + float(
                    listToStr[i_num + 1]
                )

        return catdict, biodict

    def updateDictValues(self, categories_temp: dict, biotypes_temp: dict):
        """
        Takes categories and biotypes data from dictonaries containg values
        from one file and add them into the the dictionaries previously
        parsed after matching with the dictionary keys.
        """
        if "temp" in self.categories:
            del self.categories["temp"]
        if "temp" in self.biotypes:
            del self.biotypes["temp"]
        del categories_temp["temp"]
        del biotypes_temp["temp"]

        # Combining both the dictionaries
        for i in categories_temp.keys():
            found = 0
            for j in self.categories.keys():
                if i == j:
                    found = 1
                    self.categories[j].update(categories_temp[i])
            if found == 0:
                self.categories[i] = categories_temp[i]

        for i in biotypes_temp.keys():
            found = 0
            for j in self.biotypes.keys():
                if i == j:
                    found = 1
                    self.biotypes[j].update(biotypes_temp[i])
            if found == 0:
                self.biotypes[i] = biotypes_temp[i]

    def print_alfa_charts(self, folder: str):
        """
        Takes in dictionary containing parsed data and passes them to MultiQC
        function to print the graphs.
        """
        self.add_section(
            name=f"{folder}-Categories",
            anchor="categories",
            plot=bargraph.plot(self.categories),
        )
        self.add_section(
            name=f"{folder}-BioTypes",
            anchor="biotypes",
            plot=bargraph.plot(self.biotypes),
        )
