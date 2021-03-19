from multiqc.utils import config


def execution_start():
    """Code to execute after the config files and
    command line flags have been parsedself.
    This setuptools hook is the earliest that will be able
    to use custom command line flags.
    """

    # Add to the search patterns used by modules
    if "ALFA" not in config.sp:
        config.update_dict(config.sp, {"ALFA": {"fn": "*ALFA_feature_counts.tsv"}})
    if "tin-score" not in config.sp:
        config.update_dict(config.sp, {"tin-score": {"fn": "TIN_score.tsv"}})
    if "zpca/pca" not in config.sp:
        config.update_dict(config.sp, {"zpca/pca": {"fn": "PCA.tsv"}})
    if "zpca/scree" not in config.sp:
        config.update_dict(config.sp, {"zpca/scree": {"fn": "scree.tsv"}})
