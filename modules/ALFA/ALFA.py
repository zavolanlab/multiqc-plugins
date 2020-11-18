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
import math

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    This class is instantiated in setup.py file and contains
    all the functions required by the ALFA plugin.
    """

    def __init__(self):
        """
        Initializes the parent function of MultiQC and calls
        relevant functions for the searching of files and the
        plotting of graphs.
        """
        super(MultiqcModule, self).__init__(
            name="ALFA",
            anchor="ALFA",
            href="https://github.com/biocompibens/ALFA",
            info="provides a global overview of features distribution \
                composing NGS dataset(s). Given a set of aligned reads \
                (BAM files) and an annotation file (GTF format with biotypes), \
                the tool produces plots of the raw and normalized distributions \
                of those reads among genomic categories (stop codon, 5'-UTR, \
                CDS, intergenic, etc.) and biotypes protein coding genes, miRNA\
                , tRNA, etc.). Whatever the sequencing technique, \
                whatever the organism.",
        )

        self.folders = []
        self.find_parent_folders()

        for folder in self.folders:
            self.reset()
            self.findLogs(folder)
            self.calculatePercentage()
            self.calculateEnrichment()
            self.print_alfa_charts(
                folder=folder, data=self.categories, kind="Categories"
            )
            self.print_alfa_charts(folder=folder, data=self.biotypes, kind="Biotypes")

    def find_parent_folders(self):
        """
        Finds the parent folders of log files to create a section for each folder
        in multiQC report.
        """
        for f in self.find_log_files("ALFA"):
            if os.path.basename(f["root"]) not in self.folders:
                self.folders = self.folders + [os.path.basename(f["root"])]

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

        self.categories_enrich = dict()
        self.categories_enrich["temp"] = {}

        self.biotypes_enrich = dict()
        self.biotypes_enrich["temp"] = {}

        self.categories_percent = dict()
        self.categories_percent["temp"] = {}

        self.biotypes_percent = dict()
        self.biotypes_percent["temp"] = {}

        self.categories_size = dict()
        self.categories_size["temp"] = {}

        self.biotypes_size = dict()
        self.biotypes_size["temp"] = {}

        self.categories_size_percent = dict()
        self.categories_size_percent["temp"] = {}

        self.biotypes_size_percent = dict()
        self.biotypes_size_percent["temp"] = {}

        self.categories_temp = dict()
        self.categories_temp["temp"] = {}

        self.biotypes_temp = dict()
        self.biotypes_temp["temp"] = {}

        self.categories_size_temp = dict()
        self.categories_size_temp["temp"] = {}

        self.biotypes_size_temp = dict()
        self.biotypes_size_temp["temp"] = {}

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
                    (
                        self.categories_temp,
                        self.biotypes_temp,
                        self.categories_size_temp,
                        self.biotypes_size_temp,
                    ) = self.parse_alfa_logs(f, filename)
                    self.updateDictValues(
                        self.categories_temp,
                        self.biotypes_temp,
                        self.categories_size_temp,
                        self.biotypes_size_temp,
                    )

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

        catdict_size = dict()
        catdict["temp"] = {}

        biodict_size = dict()
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

                catdict_size.setdefault(filename, {})
                biodict_size.setdefault(filename, {})

                catdict_size[filename].setdefault(cat, 0)
                biodict_size[filename].setdefault(bio, 0)

                catdict[filename][cat] = catdict[filename][cat] + float(
                    listToStr[i_num + 1]
                )
                biodict[filename][bio] = biodict[filename][bio] + float(
                    listToStr[i_num + 1]
                )

                catdict_size[filename][cat] = catdict_size[filename][cat] + float(
                    listToStr[i_num + 2]
                )
                biodict_size[filename][bio] = biodict_size[filename][bio] + float(
                    listToStr[i_num + 2]
                )
        return catdict, biodict, catdict_size, biodict_size

    def updateDictValues(
        self,
        categories_temp: dict,
        biotypes_temp: dict,
        categories_size_temp: dict,
        biotypes_size_temp: dict,
    ):
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

        if "temp" in self.categories_size:
            del self.categories_size["temp"]
        if "temp" in self.biotypes_size:
            del self.biotypes_size["temp"]
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

        for i in categories_size_temp.keys():
            found = 0
            for j in self.categories_size.keys():
                if i == j:
                    found = 1
                    self.categories_size[j].update(categories_size_temp[i])
            if found == 0:
                self.categories_size[i] = categories_size_temp[i]

        for i in biotypes_size_temp.keys():
            found = 0
            for j in self.biotypes_size.keys():
                if i == j:
                    found = 1
                    self.biotypes_size[j].update(biotypes_size_temp[i])
            if found == 0:
                self.biotypes_size[i] = biotypes_size_temp[i]

    def print_alfa_charts(self, folder: str, data: dict, kind: str):
        """
        Takes in dictionary containing parsed data and passes them to MultiQC
        function to print the graphs.
        """
        self.add_section(
            name=f"{folder}-{kind}",
            anchor=f"{kind}",
            plot=bargraph.plot(data),
        )

        self.print_alfa_charts_enrichment(folder=folder, kind=kind)

    def print_alfa_charts_enrichment(self, folder: str, kind: str):
        """
        Takes in dictionary containing parsed data and passes them to MultiQC
        function to print the enrichment graphs.
        """
        if kind == "Categories":
            data = self.categories_enrich
        else:
            data = self.biotypes_enrich

        cats = []
        for i in data.keys():
            for j in data[i]:
                cats.append(str(j))

        config = {
            # Building the plot
            "id": "<random string>",  # HTML ID used for plot
            "cpswitch": False,  # Show the 'Counts / Percentages' switch?
            "cpswitch_c_active": True,
            # Initial display with 'Counts' specified? False for percentages.
            "cpswitch_counts_label": "Counts",  # Label for 'Counts' button
            "cpswitch_percent_label": "Percentages",  # Label for 'Percentages' button
            "logswitch": False,  # Show the 'Log10' switch?
            "logswitch_active": False,  # Initial display with 'Log10' active?
            "logswitch_label": "Log10",  # Label for 'Log10' button
            "hide_zero_cats": False,  # Hide categories where data for all samples is 0
            # Customising the plot
            "title": None,  # Plot title - should be in format "Module Name: Plot Title"
            "xlab": None,  # X axis label
            "ylab": None,  # Y axis label
            "ymax": None,  # Max y limit
            "ymin": None,  # Min y limit
            "yCeiling": None,
            # Maximum value for automatic axis limit (good for percentages)
            "yFloor": None,  # Minimum value for automatic axis limit
            "yMinRange": None,  # Minimum range for axis
            "yDecimals": True,  # Set to false to only show integer labels
            "ylab_format": None,  # Format string for x axis labels. Defaults to {value}
            "stacking": None,  # Set to None to have category bars side by side
            "use_legend": True,  # Show / hide the legend
            "click_func": None,
            # Javascript function to be called when a point is clicked
            "cursor": None,  # CSS mouse cursor type.
            "tt_decimals": 2,  # Number of decimal places to use in the tooltip number
            "tt_suffix": "",  # Suffix to add after tooltip number
            "tt_percentages": False,
            # Show the percentages of each count in the tooltip
        }
        self.add_section(
            name=f"{folder}-{kind}-Enrichment",
            anchor=f"{kind}",
            plot=bargraph.plot(data, cats, config),
        )

    def calculatePercentage(self):
        """
        Calculates the percentage of each type, needed for the enrichment graphs.
        """
        if "temp" in self.categories_percent:
            del self.categories_percent["temp"]
        if "temp" in self.biotypes_percent:
            del self.biotypes_percent["temp"]
        if "temp" in self.categories_size_percent:
            del self.categories_size_percent["temp"]
        if "temp" in self.biotypes_size_percent:
            del self.biotypes_size_percent["temp"]

        categories_total = dict()
        biotypes_total = dict()
        for i in self.categories.keys():
            for j in self.categories[i]:
                categories_total[i] = categories_total.get(i, 0) + self.categories[i][j]
        for i in self.biotypes.keys():
            for j in self.biotypes[i]:
                biotypes_total[i] = biotypes_total.get(i, 0) + self.biotypes[i][j]

        # Divinding by total subiotypes
        for i in self.categories.keys():
            self.categories_percent.setdefault(i, {})
            for j in self.categories[i]:
                self.categories_percent[i].setdefault(j, 0)
                self.categories_percent[i][j] = (
                    self.categories[i][j] / categories_total[i] * 100
                )
        for i in self.biotypes.keys():
            self.biotypes_percent.setdefault(i, {})
            for j in self.biotypes[i]:
                self.biotypes_percent[i].setdefault(j, 0)
                self.biotypes_percent[i][j] = (
                    self.biotypes[i][j] / biotypes_total[i] * 100
                )

        categories_size_total = dict()
        biotypes_size_total = dict()
        for i in self.categories_size.keys():
            for j in self.categories_size[i]:
                categories_size_total[i] = (
                    categories_size_total.get(i, 0) + self.categories_size[i][j]
                )
        for i in self.biotypes_size.keys():
            for j in self.biotypes_size[i]:
                biotypes_size_total[i] = (
                    biotypes_size_total.get(i, 0) + self.biotypes_size[i][j]
                )

        # Divinding by total subiotypes
        for i in self.categories_size.keys():
            self.categories_size_percent.setdefault(i, {})
            for j in self.categories_size[i]:
                self.categories_size_percent[i].setdefault(j, 0)
                self.categories_size_percent[i][j] = (
                    self.categories_size[i][j] / categories_size_total[i] * 100
                )
        for i in self.biotypes_size.keys():
            self.biotypes_size_percent.setdefault(i, {})
            for j in self.biotypes_size[i]:
                self.biotypes_size_percent[i].setdefault(j, 0)
                self.biotypes_size_percent[i][j] = (
                    self.biotypes_size[i][j] / biotypes_size_total[i] * 100
                )

    def calculateEnrichment(self):
        """
        Calculates percentage of nucleotides divided by percentage of genome for
        total enrichment or depletion.
        """
        if "temp" in self.categories_enrich:
            del self.categories_enrich["temp"]
        if "temp" in self.biotypes_enrich:
            del self.biotypes_enrich["temp"]

        for i in self.categories.keys():
            self.categories_enrich.setdefault(i, {})
            for j in self.categories[i]:
                self.categories_enrich[i].setdefault(j, 0)
                self.categories_enrich[i][j] = math.log2(
                    self.categories_percent[i][j] / self.categories_size_percent[i][j]
                )

        for i in self.biotypes.keys():
            self.biotypes_enrich.setdefault(i, {})
            for j in self.biotypes[i]:
                self.biotypes_enrich[i].setdefault(j, 0)
                self.biotypes_enrich[i][j] = math.log2(
                    self.biotypes_percent[i][j] / self.biotypes_size_percent[i][j]
                )
