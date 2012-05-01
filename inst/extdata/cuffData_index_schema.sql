CREATE INDEX "genes.gsn_index" ON "genes"("gene_short_name");
CREATE INDEX "genes.cc_index" ON "genes"("class_code");
CREATE INDEX "TSS.fk_TSS_genes1" ON "TSS"("gene_id");
CREATE INDEX "TSSData.fk_TSSData_TSS1" ON "TSSData"("TSS_group_id");
CREATE INDEX "TSSData.fk_TSSData_samples1" ON "TSSData"("sample_name");
CREATE INDEX "CDS.fk_CDS_genes1" ON "CDS"("gene_id");
CREATE INDEX "CDS.fk_CDS_TSS1" ON "CDS"("TSS_group_id");
CREATE INDEX "CDSData.fk_CDSData_CDS1" ON "CDSData"("CDS_id");
CREATE INDEX "CDSData.fk_CDSData_samples1" ON "CDSData"("sample_name");
CREATE INDEX "splicingDiffData.fk_splicingDiffData_samples1" ON "splicingDiffData"("sample_1");
CREATE INDEX "splicingDiffData.fk_splicingDiffData_samples2" ON "splicingDiffData"("sample_2");
CREATE INDEX "splicingDiffData.fk_splicingDiffData_TSS1" ON "splicingDiffData"("TSS_group_id");
CREATE INDEX "splicingDiffData.fk_splicingDiffData_genes1" ON "splicingDiffData"("gene_id");
CREATE INDEX "TSSExpDiffData.fk_TSSExpDiffData_TSS1" ON "TSSExpDiffData"("TSS_group_id");
CREATE INDEX "TSSExpDiffData.fk_TSSExpDiffData_samples1" ON "TSSExpDiffData"("sample_1");
CREATE INDEX "TSSExpDiffData.fk_TSSExpDiffData_samples2" ON "TSSExpDiffData"("sample_2");
CREATE INDEX "TSSExpDiffData.TSSExpDiffData_sig_index" ON "TSSExpDiffData"("test_stat","p_value","q_value","significant");
CREATE INDEX "CDSDiffData.fk_CDSDiffData_samples1" ON "CDSDiffData"("sample_1");
CREATE INDEX "CDSDiffData.fk_CDSDiffData_samples2" ON "CDSDiffData"("sample_2");
CREATE INDEX "CDSDiffData.fk_CDSDiffData_genes1" ON "CDSDiffData"("gene_id");
CREATE INDEX "CDSExpDiffData.fk_CDSExpDiffData_CDS1" ON "CDSExpDiffData"("CDS_id");
CREATE INDEX "CDSExpDiffData.fk_CDSExpDiffData_samples1" ON "CDSExpDiffData"("sample_1");
CREATE INDEX "CDSExpDiffData.fk_CDSExpDiffData_samples2" ON "CDSExpDiffData"("sample_2");
CREATE INDEX "CDSExpDiffData.CDSExpDiffData_sig_index" ON "CDSExpDiffData"("test_stat","p_value","q_value","significant");
CREATE INDEX "promoterDiffData.fk_promoterDiffData_genes1" ON "promoterDiffData"("gene_id");
CREATE INDEX "promoterDiffData.fk_promoterDiffData_samples1" ON "promoterDiffData"("sample_1");
CREATE INDEX "promoterDiffData.fk_promoterDiffData_samples2" ON "promoterDiffData"("sample_2");
CREATE INDEX "geneFeatures.fk_geneFeatures_genes1" ON "geneFeatures"("gene_id");
CREATE INDEX "TSSFeatures.fk_TSSFeatures_TSS1" ON "TSSFeatures"("TSS_group_id");
CREATE INDEX "CDSFeatures.fk_CDSFeatures_CDS1" ON "CDSFeatures"("CDS_id");
CREATE INDEX "geneCount.fk_geneCount_samples1" ON "geneCount"("sample_name");
CREATE INDEX "geneCount.fk_geneCount_genes1" ON "geneCount"("gene_id");
CREATE INDEX "CDSCount.fk_CDSCount_CDS1" ON "CDSCount"("CDS_id");
CREATE INDEX "CDSCount.fk_CDSCount_samples1" ON "CDSCount"("sample_name");
CREATE INDEX "TSSCount.fk_TSSCount_TSS1" ON "TSSCount"("TSS_group_id");
CREATE INDEX "TSSCount.fk_TSSCount_samples1" ON "TSSCount"("sample_name");
CREATE INDEX "replicates.fk_replicates_samples1" ON "replicates"("sample_name");
CREATE INDEX "geneReplicateData.fk_geneData_genes1" ON "geneReplicateData"("gene_id");
CREATE INDEX "geneReplicateData.fk_geneReplicateData_replicates1" ON "geneReplicateData"("rep_name");
CREATE INDEX "CDSReplicateData.fk_CDSReplicateData_replicates1" ON "CDSReplicateData"("rep_name");
CREATE INDEX "CDSReplicateData.fk_CDSReplicateData_CDS1" ON "CDSReplicateData"("CDS_id");
CREATE INDEX "TSSReplicateData.fk_TSSReplicateData_replicates1" ON "TSSReplicateData"("rep_name");
CREATE INDEX "TSSReplicateData.fk_TSSReplicateData_TSS1" ON "TSSReplicateData"("TSS_group_id");
CREATE INDEX "geneData.fk_geneData_genes1" ON "geneData"("gene_id");
CREATE INDEX "geneData.fk_geneData_samples1" ON "geneData"("sample_name");
CREATE INDEX "phenoData.fk_phenoData_samples" ON "phenoData"("sample_name");
CREATE INDEX "geneExpDiffData.fk_geneExpDiffData_genes1" ON "geneExpDiffData"("gene_id");
CREATE INDEX "geneExpDiffData.fk_geneExpDiffData_samples1" ON "geneExpDiffData"("sample_1");
CREATE INDEX "geneExpDiffData.fk_geneExpDiffData_samples2" ON "geneExpDiffData"("sample_2");
CREATE INDEX "geneExpDiffData.geneExpDiff_status_index" ON "geneExpDiffData"("status");
CREATE INDEX "geneExpDiffData.geneExpDiff_sig_index" ON "geneExpDiffData"("significant","p_value","q_value","test_stat");
CREATE INDEX "isoforms.fk_isoforms_TSS1" ON "isoforms"("TSS_group_id");
CREATE INDEX "isoforms.fk_isoforms_CDS1" ON "isoforms"("CDS_id");
CREATE INDEX "isoforms.fk_isoforms_genes1" ON "isoforms"("gene_id");
CREATE INDEX "isoformData.fk_isoformData_samples1" ON "isoformData"("sample_name");
CREATE INDEX "isoformData.fk_isoformData_isoforms1" ON "isoformData"("isoform_id");
CREATE INDEX "isoformExpDiffData.fk_isoformExpDiffData_isoforms1" ON "isoformExpDiffData"("isoform_id");
CREATE INDEX "isoformExpDiffData.fk_isoformExpDiffData_samples1" ON "isoformExpDiffData"("sample_1");
CREATE INDEX "isoformExpDiffData.fk_isoformExpDiffData_samples2" ON "isoformExpDiffData"("sample_2");
CREATE INDEX "isoformExpDiffData.isoformExpDiffData_sig_index" ON "isoformExpDiffData"("test_stat","p_value","q_value","significant");
CREATE INDEX "isoformFeatures.fk_isoformFeatures_isoforms1" ON "isoformFeatures"("isoform_id");
CREATE INDEX "features.features_seqname_index" ON "features"("seqname");
CREATE INDEX "features.features_type_index" ON "features"("type_id");
CREATE INDEX "features.features_strand_index" ON "features"("strand");
CREATE INDEX "features.features_start_end_index" ON "features"("start","end");
CREATE INDEX "features.fk_features_genes1" ON "features"("genes_gene_id");
CREATE INDEX "features.fk_features_isoforms1" ON "features"("isoforms_isoform_id");
CREATE INDEX "attributes.fk_attributes_feature_id" ON "attributes"("feature_id");
CREATE INDEX "attributes.attributes_attribute_index" ON "attributes"("attribute");
CREATE INDEX "attributes.attributes_value_index" ON "attributes"("value");
CREATE INDEX "isoformCount.fk_isoformCount_isoforms1" ON "isoformCount"("isoform_id");
CREATE INDEX "isoformCount.fk_isoformCount_samples1" ON "isoformCount"("sample_name");
CREATE INDEX "isoformReplicateData.fk_isoformReplicateData_replicates1" ON "isoformReplicateData"("rep_name");
CREATE INDEX "isoformReplicateData.fk_isoformReplicateData_isoforms1" ON "isoformReplicateData"("isoform_id");
