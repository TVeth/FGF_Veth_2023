library(artMS)

artmsQualityControlEvidenceExtended(
  evidence_file = "evidence.txt",
  keys_file = "keys.txt")

artmsQualityControlSummaryExtended(summary_file = "summary.txt",
                                   keys_file = "keys.txt")

artmsQuantification(yaml_config_file = "my_config.yaml")

artmsAnalysisQuantifications(log2fc_file = "results.txt",
                             modelqc_file = "results_ModelQC.txt",
                             species = "human",
                             output_dir = "AnalysisQuantifications")