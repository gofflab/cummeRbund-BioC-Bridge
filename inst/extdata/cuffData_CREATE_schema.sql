-- Creator:       MySQL Workbench 5.2.33/ExportSQLite plugin 2009.12.02
-- Author:        Loyal Goff
-- Caption:       New Model
-- Project:       Name of the project
-- Changed:       2012-04-30 22:21
-- Created:       2011-05-02 12:52
PRAGMA foreign_keys = OFF;

-- Schema: cuffData
BEGIN;
DROP TABLE IF EXISTS "genes";
CREATE TABLE "genes"(
  "gene_id" VARCHAR(45) PRIMARY KEY NOT NULL,
  "class_code" VARCHAR(45),
  "nearest_ref_id" VARCHAR(45),
  "gene_short_name" VARCHAR(45),
  "locus" VARCHAR(45),
  "length" INTEGER,
  "coverage" FLOAT
);
DROP TABLE IF EXISTS "biasData";
CREATE TABLE "biasData"(
  "biasData_id" INTEGER PRIMARY KEY NOT NULL
);
DROP TABLE IF EXISTS "samples";
CREATE TABLE "samples"(
  "sample_index" INTEGER NOT NULL,
  "sample_name" VARCHAR(45) PRIMARY KEY NOT NULL
);
DROP TABLE IF EXISTS "TSS";
CREATE TABLE "TSS"(
  "TSS_group_id" VARCHAR(45) PRIMARY KEY NOT NULL,
  "class_code" VARCHAR(45),
  "nearest_ref_id" VARCHAR(45),
  "gene_id" VARCHAR(45) NOT NULL,
  "locus" VARCHAR(45),
  "length" INTEGER,
  "coverage" FLOAT,
  CONSTRAINT "fk_TSS_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
DROP TABLE IF EXISTS "TSSData";
CREATE TABLE "TSSData"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "fpkm" FLOAT,
  "conf_hi" FLOAT,
  "conf_lo" FLOAT,
  "quant_status" VARCHAR(45),
  CONSTRAINT "fk_TSSData_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id"),
  CONSTRAINT "fk_TSSData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "CDS";
CREATE TABLE "CDS"(
  "CDS_id" VARCHAR(45) PRIMARY KEY NOT NULL,
  "class_code" VARCHAR(45),
  "nearest_ref_id" VARCHAR(45),
  "gene_id" VARCHAR(45),
  "TSS_group_id" VARCHAR(45),
  "locus" VARCHAR(45),
  "length" INTEGER,
  "coverage" FLOAT,
  CONSTRAINT "fk_CDS_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_CDS_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id")
);
DROP TABLE IF EXISTS "CDSData";
CREATE TABLE "CDSData"(
  "CDS_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "fpkm" FLOAT,
  "conf_hi" FLOAT,
  "conf_lo" FLOAT,
  "quant_status" VARCHAR(45),
  CONSTRAINT "fk_CDSData_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id"),
  CONSTRAINT "fk_CDSData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "splicingDiffData";
CREATE TABLE "splicingDiffData"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "JS_dist" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_splicingDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_splicingDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_splicingDiffData_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id"),
  CONSTRAINT "fk_splicingDiffData_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
DROP TABLE IF EXISTS "TSSExpDiffData";
CREATE TABLE "TSSExpDiffData"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "ln_fold_change" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_TSSExpDiffData_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id"),
  CONSTRAINT "fk_TSSExpDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_TSSExpDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "CDSDiffData";
CREATE TABLE "CDSDiffData"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "JS_dist" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_CDSDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_CDSDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_CDSDiffData_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
DROP TABLE IF EXISTS "CDSExpDiffData";
CREATE TABLE "CDSExpDiffData"(
  "CDS_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "ln_fold_change" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_CDSExpDiffData_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id"),
  CONSTRAINT "fk_CDSExpDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_CDSExpDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "promoterDiffData";
CREATE TABLE "promoterDiffData"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "JS_dist" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_promoterDiffData_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_promoterDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_promoterDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "geneFeatures";
CREATE TABLE "geneFeatures"(
  "gene_id" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_geneFeatures_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
DROP TABLE IF EXISTS "TSSFeatures";
CREATE TABLE "TSSFeatures"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_TSSFeatures_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id")
);
DROP TABLE IF EXISTS "CDSFeatures";
CREATE TABLE "CDSFeatures"(
  "CDS_id" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_CDSFeatures_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id")
);
DROP TABLE IF EXISTS "model_transcripts";
CREATE TABLE "model_transcripts"(
  "model_transcript_id" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL
);
DROP TABLE IF EXISTS "geneCount";
CREATE TABLE "geneCount"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "count" FLOAT,
  "variance" FLOAT,
  "status" VARCHAR(45),
  CONSTRAINT "fk_geneCount_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_geneCount_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
DROP TABLE IF EXISTS "CDSCount";
CREATE TABLE "CDSCount"(
  "CDS_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "count" FLOAT,
  "variance" FLOAT,
  "status" VARCHAR(45),
  CONSTRAINT "fk_CDSCount_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id"),
  CONSTRAINT "fk_CDSCount_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "TSSCount";
CREATE TABLE "TSSCount"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "count" FLOAT,
  "variance" FLOAT,
  "status" VARCHAR(45),
  CONSTRAINT "fk_TSSCount_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id"),
  CONSTRAINT "fk_TSSCount_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "replicates";
CREATE TABLE "replicates"(
  "file" INTEGER NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "rep_name" VARCHAR(45) PRIMARY KEY NOT NULL,
  "total_mass" FLOAT,
  "norm_mass" FLOAT,
  "internal_scale" FLOAT,
  "external_scale" FLOAT,
  CONSTRAINT "fk_replicates_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "geneReplicateData";
CREATE TABLE "geneReplicateData"(
  "gene_id" VARCHAR(45) NOT NULL,
  "rep_name" VARCHAR(45) NOT NULL,
  "raw_frags" FLOAT,
  "internal_scaled_frags" FLOAT,
  "external_scaled_frags" FLOAT,
  "fpkm" FLOAT,
  "effective_length" FLOAT,
  "status" VARCHAR(45),
  CONSTRAINT "fk_geneData_genes10"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_geneReplicateData_replicates1"
    FOREIGN KEY("rep_name")
    REFERENCES "replicates"("rep_name")
);
DROP TABLE IF EXISTS "CDSReplicateData";
CREATE TABLE "CDSReplicateData"(
  "CDS_id" VARCHAR(45) NOT NULL,
  "rep_name" VARCHAR(45) NOT NULL,
  "raw_frags" FLOAT,
  "internal_scaled_frags" FLOAT,
  "external_scaled_frags" FLOAT,
  "fpkm" FLOAT,
  "effective_length" FLOAT,
  "status" VARCHAR(45),
  CONSTRAINT "fk_geneReplicateData_replicates100"
    FOREIGN KEY("rep_name")
    REFERENCES "replicates"("rep_name"),
  CONSTRAINT "fk_CDSReplicateData_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id")
);
DROP TABLE IF EXISTS "TSSReplicateData";
CREATE TABLE "TSSReplicateData"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  "rep_name" VARCHAR(45) NOT NULL,
  "raw_frags" FLOAT,
  "internal_scaled_frags" FLOAT,
  "external_scaled_frags" FLOAT,
  "fpkm" FLOAT,
  "effective_length" FLOAT,
  "status" VARCHAR(45),
  CONSTRAINT "fk_geneReplicateData_replicates10000"
    FOREIGN KEY("rep_name")
    REFERENCES "replicates"("rep_name"),
  CONSTRAINT "fk_TSSReplicateData_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id")
);
DROP TABLE IF EXISTS "runInfo";
CREATE TABLE "runInfo"(
  "runInfo_id" INTEGER PRIMARY KEY NOT NULL,
  "param" VARCHAR(45),
  "value" TEXT
);
DROP TABLE IF EXISTS "geneData";
CREATE TABLE "geneData"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "fpkm" FLOAT,
  "conf_hi" FLOAT,
  "conf_lo" FLOAT,
  "quant_status" VARCHAR(45),
  CONSTRAINT "fk_geneData_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_geneData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "phenoData";
CREATE TABLE "phenoData"(
  "sample_name" VARCHAR(45) NOT NULL,
  "parameter" VARCHAR(45) NOT NULL,
  "value" VARCHAR(45),
  CONSTRAINT "fk_phenoData_samples"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "geneExpDiffData";
CREATE TABLE "geneExpDiffData"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "ln_fold_change" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_geneExpDiffData_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_geneExpDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_geneExpDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "isoforms";
CREATE TABLE "isoforms"(
  "isoform_id" VARCHAR(45) PRIMARY KEY NOT NULL,
  "gene_id" VARCHAR(45),
  "CDS_id" VARCHAR(45),
  "TSS_group_id" VARCHAR(45),
  "class_code" VARCHAR(45),
  "nearest_ref_id" VARCHAR(45),
  "locus" VARCHAR(45),
  "length" INTEGER,
  "coverage" FLOAT,
  CONSTRAINT "fk_isoforms_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id"),
  CONSTRAINT "fk_isoforms_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id"),
  CONSTRAINT "fk_isoforms_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
DROP TABLE IF EXISTS "isoformData";
CREATE TABLE "isoformData"(
  "isoform_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "fpkm" FLOAT NOT NULL,
  "conf_hi" FLOAT,
  "conf_lo" FLOAT,
  "quant_status" VARCHAR(45),
  CONSTRAINT "fk_isoformData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_isoformData_isoforms1"
    FOREIGN KEY("isoform_id")
    REFERENCES "isoforms"("isoform_id")
);
DROP TABLE IF EXISTS "isoformExpDiffData";
CREATE TABLE "isoformExpDiffData"(
  "isoform_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "ln_fold_change" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_isoformExpDiffData_isoforms1"
    FOREIGN KEY("isoform_id")
    REFERENCES "isoforms"("isoform_id"),
  CONSTRAINT "fk_isoformExpDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_isoformExpDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "isoformFeatures";
CREATE TABLE "isoformFeatures"(
  "isoform_id" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_isoformFeatures_isoforms1"
    FOREIGN KEY("isoform_id")
    REFERENCES "isoforms"("isoform_id")
);
DROP TABLE IF EXISTS "features";
CREATE TABLE "features"(
--   GTF Features (all lines/records from reference .gtf file)
  "feature_id" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  "genes_gene_id" VARCHAR(45) NOT NULL,
  "isoforms_isoform_id" VARCHAR(45) NOT NULL,
  "seqname" VARCHAR(45) NOT NULL,
  "source" VARCHAR(45) NOT NULL,
  "type_id" INTEGER,
  "start" INTEGER,
  "end" INTEGER,
  "score" FLOAT,
  "strand" VARCHAR(45),
  "frame" VARCHAR(45),
  CONSTRAINT "fk_features_genes1"
    FOREIGN KEY("genes_gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_features_isoforms1"
    FOREIGN KEY("isoforms_isoform_id")
    REFERENCES "isoforms"("isoform_id")
);
DROP TABLE IF EXISTS "attribtes";
CREATE TABLE "attributes"(
  "attribute_lookup_id" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  "feature_id" INTEGER NOT NULL,
  "attribute" VARCHAR(45) NOT NULL,
  "value" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_attribute_lookup_features1"
    FOREIGN KEY("feature_id")
    REFERENCES "features"("feature_id")
);
DROP TABLE IF EXISTS "isoformCount";
CREATE TABLE "isoformCount"(
  "isoform_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "count" FLOAT,
  "variance" FLOAT,
  "status" VARCHAR(45),
  CONSTRAINT "fk_isoformCount_isoforms1"
    FOREIGN KEY("isoform_id")
    REFERENCES "isoforms"("isoform_id"),
  CONSTRAINT "fk_isoformCount_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "isoformReplicateData";
CREATE TABLE "isoformReplicateData"(
  "isoform_id" VARCHAR(45) NOT NULL,
  "rep_name" VARCHAR(45) NOT NULL,
  "raw_frags" FLOAT,
  "internal_scaled_frags" FLOAT,
  "external_scaled_frags" FLOAT,
  "fpkm" FLOAT,
  "effective_length" FLOAT,
  "status" VARCHAR(45),
  CONSTRAINT "fk_geneReplicateData_replicates10"
    FOREIGN KEY("rep_name")
    REFERENCES "replicates"("rep_name"),
  CONSTRAINT "fk_isoformReplicateData_isoforms1"
    FOREIGN KEY("isoform_id")
    REFERENCES "isoforms"("isoform_id")
);
COMMIT;
