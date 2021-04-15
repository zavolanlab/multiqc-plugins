###############################################################################
#
#   Main script which does the job of parsing values and plotting graphs
#   for zpca values
#
#   AUTHOR: Krish Agarwal
#   AFFILIATION: University_of_Basel
#   CONTACT: akrish136@gmail.com
#   CREATED: 15-03-2021
#   LICENSE: Apache_2.0
#
###############################################################################

from __future__ import print_function
import logging

from multiqc.plots import scatter
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    This class is instantiated in setup.py file and contains
    all the functions required by the zpca plugin.
    """

    def __init__(self):
        """
        Initializes the parent function of MultiQC and calls
        relevant functions for the searching of files and the
        plotting of graphs.
        """

        super(MultiqcModule, self).__init__(
            name="zpca",
            anchor="zpca",
            href="https://github.com/zavolanlab/zpca",
            info="- PCA analysis",
        )

        self.number = 0
        self.findLogs()

    def findLogs(self):
        """
        Find files matching with the regex in and passes them to the
        parsing function one by one and also plots the graph in the end.
        """
        for f in self.find_log_files("zpca/scree"):
            self.number += 1
            data_scree = self.parse_scree_logs(f)

            exp_car_str = "Percentage of Explained Variance"
            for f1 in self.find_log_files("zpca/pca"):
                if f["root"] == f1["root"]:
                    self.number += 1
                    data_list = self.parse_zpca_logs(f1)
                    if len(data_list) == 3:
                        data_pc1_pc2, data_pc1_pc3, data_pc2_pc3 = (
                            data_list[0],
                            data_list[1],
                            data_list[2],
                        )
                        config_pc1_pc2 = {
                            "xlab": (
                                f"PC1 ({data_scree['PC1'][exp_car_str]}% variance"
                                " explained)"
                            ),
                            "ylab": (
                                f"PC2 ({data_scree['PC2'][exp_car_str]}% variance"
                                " explained)"
                            ),
                        }
                        config_pc1_pc3 = {
                            "xlab": (
                                f"PC1 ({data_scree['PC1'][exp_car_str]}% variance"
                                " explained)"
                            ),
                            "ylab": (
                                f"PC3 ({data_scree['PC3'][exp_car_str]}% variance"
                                " explained)"
                            ),
                        }
                        config_pc2_pc3 = {
                            "xlab": (
                                f"PC2 ({data_scree['PC2'][exp_car_str]}% variance"
                                " explained)"
                            ),
                            "ylab": (
                                f"PC3 ({data_scree['PC3'][exp_car_str]}% variance"
                                " explained)"
                            ),
                        }

                        self.add_section(
                            name="PCA components: 1 & 2",
                            anchor="zpca",
                            plot=scatter.plot(data_pc1_pc2, config_pc1_pc2),
                        )

                        self.add_section(
                            name="PCA components: 1 & 3",
                            anchor="zpca",
                            plot=scatter.plot(data_pc1_pc3, config_pc1_pc3),
                        )

                        self.add_section(
                            name="PCA components: 2 & 3",
                            anchor="zpca",
                            plot=scatter.plot(data_pc2_pc3, config_pc2_pc3),
                        )
                    elif len(data_list) == 1:
                        data_pc1_pc2 = data_list[0]
                        config_pc1_pc2 = {
                            "xlab": (
                                f"PC1 ({data_scree['PC1'][exp_car_str]}% variance"
                                " explained)"
                            ),
                            "ylab": (
                                f"PC2 ({data_scree['PC2'][exp_car_str]}% variance"
                                " explained)"
                            ),
                        }

                        self.add_section(
                            name="PCA components: 1 & 2",
                            anchor="zpca",
                            plot=scatter.plot(data_pc1_pc2, config_pc1_pc2),
                        )

        if self.number == 0:
            raise UserWarning

    def parse_scree_logs(self, f: dict):
        """
        Parses Scree.tsv and returns explained variance
        in percentage of principal components in a dictionary.
        """
        word = []
        words = []
        listToStr = []

        # Using \t as a seperator for elements in tsv file
        for char in f["f"]:
            if char != "\t" and char != "\n":
                word.append(char)
            else:
                words.append(word)
                word = []

        # Adding each element of the tsv file in a list
        for k in words:
            listToStr.append("".join([str(elem) for elem in k]))

        """
        Appropriate function is called depending upon whether
        the required components (PC_3_comp_list, PC_2_comp_list)
        are present in tsv file(listToStr).
        """
        PC_3_comp_list = ["PC1", "PC2", "PC3"]
        PC_2_comp_list = ["PC1", "PC2"]
        if all(elem in listToStr for elem in PC_3_comp_list):
            return self._parse_scree_helper_3(listToStr)
        if all(elem in listToStr for elem in PC_2_comp_list):
            return self._parse_scree_helper_2(listToStr)
        else:
            return {}

    def _parse_scree_helper_3(self, listToStr):
        """
        Helper function which will parse scree logs for 3 components.
        """
        exp_car_str = "Percentage of Explained Variance"
        scree_data = {
            "PC1": {exp_car_str: float(listToStr[-3])},
            "PC2": {exp_car_str: float(listToStr[-2])},
            "PC3": {exp_car_str: float(listToStr[-1])},
        }
        return scree_data

    def _parse_scree_helper_2(self, listToStr):
        """
        Helper function which will parse scree logs for 2 components.
        """
        exp_car_str = "Percentage of Explained Variance"
        scree_data = {
            "PC1": {exp_car_str: float(listToStr[-2])},
            "PC2": {exp_car_str: float(listToStr[-1])},
        }
        return scree_data

    def parse_zpca_logs(self, f: dict):
        """
        Parses the PCA.tsv file by storing the values of ZPCA in
        a dictionary and returns three dictionaries which will
        directly be used to plot the graphs in findLogs().
        """
        word = []
        words = []
        listToStr = []

        # Using \t as a seperator for elements in tsv file
        for char in f["f"]:
            if char != "\t" and char != "\n":
                word.append(char)
            else:
                words.append(word)
                word = []

        # Adding each element of the tsv file in a list
        for k in words:
            listToStr.append("".join([str(elem) for elem in k]))

        """
        Appropriate function is called depending upon whether
        the required components (PC_3_comp_list, PC_2_comp_list)
        are present in tsv file(listToStr).
        """
        PC_3_comp_list = ["PC1", "PC2", "PC3"]
        PC_2_comp_list = ["PC1", "PC2"]
        if all(elem in listToStr for elem in PC_3_comp_list):
            return self._parse_zpca_helper_3(listToStr)
        if all(elem in listToStr for elem in PC_2_comp_list):
            return self._parse_zpca_helper_2(listToStr)
        else:
            return []

    def _parse_zpca_helper_3(self, listToStr: list):
        """
        Helper function which will parse logs with 3 components.
        """
        data_dict = {}

        for i_num in range(0, len(listToStr), 4):
            if i_num == 0:
                continue
            data_dict[listToStr[i_num]] = [
                listToStr[i_num + 1],
                listToStr[i_num + 2],
                listToStr[i_num + 3],
            ]

        data_pc1_pc2 = dict()
        data_pc1_pc3 = dict()
        data_pc2_pc3 = dict()

        for i in data_dict:
            data_pc1_pc2[i] = {
                "x": float(data_dict[i][0]),
                "y": float(data_dict[i][1]),
                "color": "#58a0c3",
            }
            data_pc1_pc3[i] = {
                "x": float(data_dict[i][0]),
                "y": float(data_dict[i][2]),
                "color": "#58a0c3",
            }
            data_pc2_pc3[i] = {
                "x": float(data_dict[i][1]),
                "y": float(data_dict[i][2]),
                "color": "#58a0c3",
            }

        return [data_pc1_pc2, data_pc1_pc3, data_pc2_pc3]

    def _parse_zpca_helper_2(self, listToStr: list):
        """
        Helper function which will parse logs with 2 components.
        """
        data_dict = {}

        for i_num in range(0, len(listToStr), 3):
            if i_num == 0:
                continue
            data_dict[listToStr[i_num]] = [
                listToStr[i_num + 1],
                listToStr[i_num + 2],
            ]

        data_pc1_pc2 = dict()

        for i in data_dict:
            data_pc1_pc2[i] = {
                "x": float(data_dict[i][0]),
                "y": float(data_dict[i][1]),
                "color": "#58a0c3",
            }

        return [data_pc1_pc2]
