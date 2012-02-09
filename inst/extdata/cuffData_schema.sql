-- Creator:       MySQL Workbench 5.2.33/ExportSQLite plugin 2009.12.02
-- Author:        Loyal Goff
-- Caption:       CuffData.db Model
-- Project:       cummeRbund
-- Changed:       2011-08-02 14:03
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
  "sample_index" INTEGER PRIMARY KEY NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL
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
CREATE INDEX "TSS.fk_TSS_genes1" ON "TSS"("gene_id");
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
CREATE INDEX "TSSData.fk_TSSData_TSS1" ON "TSSData"("TSS_group_id");
CREATE INDEX "TSSData.fk_TSSData_samples1" ON "TSSData"("sample_name");
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
CREATE INDEX "CDS.fk_CDS_genes1" ON "CDS"("gene_id");
CREATE INDEX "CDS.fk_CDS_TSS1" ON "CDS"("TSS_group_id");
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
CREATE INDEX "CDSData.fk_CDSData_CDS1" ON "CDSData"("CDS_id");
CREATE INDEX "CDSData.fk_CDSData_samples1" ON "CDSData"("sample_name");
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
CREATE INDEX "splicingDiffData.fk_splicingDiffData_samples1" ON "splicingDiffData"("sample_1");
CREATE INDEX "splicingDiffData.fk_splicingDiffData_samples2" ON "splicingDiffData"("sample_2");
CREATE INDEX "splicingDiffData.fk_splicingDiffData_TSS1" ON "splicingDiffData"("TSS_group_id");
CREATE INDEX "splicingDiffData.fk_splicingDiffData_genes1" ON "splicingDiffData"("gene_id");
DROP TABLE IF EXISTS "TSSExpDiffData";
CREATE TABLE "TSSExpDiffData"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "log2_fold_change" FLOAT,
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
CREATE INDEX "TSSExpDiffData.fk_TSSExpDiffData_TSS1" ON "TSSExpDiffData"("TSS_group_id");
CREATE INDEX "TSSExpDiffData.fk_TSSExpDiffData_samples1" ON "TSSExpDiffData"("sample_1");
CREATE INDEX "TSSExpDiffData.fk_TSSExpDiffData_samples2" ON "TSSExpDiffData"("sample_2");
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
CREATE INDEX "CDSDiffData.fk_CDSDiffData_samples1" ON "CDSDiffData"("sample_1");
CREATE INDEX "CDSDiffData.fk_CDSDiffData_samples2" ON "CDSDiffData"("sample_2");
CREATE INDEX "CDSDiffData.fk_CDSDiffData_genes1" ON "CDSDiffData"("gene_id");
DROP TABLE IF EXISTS "CDSExpDiffData";
CREATE TABLE "CDSExpDiffData"(
  "CDS_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "log2_fold_change" FLOAT,
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
CREATE INDEX "CDSExpDiffData.fk_CDSExpDiffData_CDS1" ON "CDSExpDiffData"("CDS_id");
CREATE INDEX "CDSExpDiffData.fk_CDSExpDiffData_samples1" ON "CDSExpDiffData"("sample_1");
CREATE INDEX "CDSExpDiffData.fk_CDSExpDiffData_samples2" ON "CDSExpDiffData"("sample_2");
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
CREATE INDEX "promoterDiffData.fk_promoterDiffData_genes1" ON "promoterDiffData"("gene_id");
CREATE INDEX "promoterDiffData.fk_promoterDiffData_samples1" ON "promoterDiffData"("sample_1");
CREATE INDEX "promoterDiffData.fk_promoterDiffData_samples2" ON "promoterDiffData"("sample_2");
DROP TABLE IF EXISTS "geneFeatures";
CREATE TABLE "geneFeatures"(
  "gene_id" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_geneFeatures_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
CREATE INDEX "geneFeatures.fk_geneFeatures_genes1" ON "geneFeatures"("gene_id");
DROP TABLE IF EXISTS "TSSFeatures";
CREATE TABLE "TSSFeatures"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_TSSFeatures_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id")
);
CREATE INDEX "TSSFeatures.fk_TSSFeatures_TSS1" ON "TSSFeatures"("TSS_group_id");
DROP TABLE IF EXISTS "CDSFeatures";
CREATE TABLE "CDSFeatures"(
  "CDS_id" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_CDSFeatures_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id")
);
CREATE INDEX "CDSFeatures.fk_CDSFeatures_CDS1" ON "CDSFeatures"("CDS_id");
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
CREATE INDEX "geneData.fk_geneData_genes1" ON "geneData"("gene_id");
CREATE INDEX "geneData.fk_geneData_samples1" ON "geneData"("sample_name");
DROP TABLE IF EXISTS "phenoData";
CREATE TABLE "phenoData"(
  "sample_name" VARCHAR(45) NOT NULL,
  "parameter" VARCHAR(45) NOT NULL,
  "value" VARCHAR(45),
  CONSTRAINT "fk_phenoData_samples"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "phenoData.fk_phenoData_samples" ON "phenoData"("sample_name");
DROP TABLE IF EXISTS "geneExpDiffData";
CREATE TABLE "geneExpDiffData"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "log2_fold_change" FLOAT,
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
CREATE INDEX "geneExpDiffData.fk_geneExpDiffData_genes1" ON "geneExpDiffData"("gene_id");
CREATE INDEX "geneExpDiffData.fk_geneExpDiffData_samples1" ON "geneExpDiffData"("sample_1");
CREATE INDEX "geneExpDiffData.fk_geneExpDiffData_samples2" ON "geneExpDiffData"("sample_2");
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
CREATE INDEX "isoforms.fk_isoforms_TSS1" ON "isoforms"("TSS_group_id");
CREATE INDEX "isoforms.fk_isoforms_CDS1" ON "isoforms"("CDS_id");
CREATE INDEX "isoforms.fk_isoforms_genes1" ON "isoforms"("gene_id");
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
CREATE INDEX "isoformData.fk_isoformData_samples1" ON "isoformData"("sample_name");
CREATE INDEX "isoformData.fk_isoformData_isoforms1" ON "isoformData"("isoform_id");
DROP TABLE IF EXISTS "isoformExpDiffData";
CREATE TABLE "isoformExpDiffData"(
  "isoform_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "log2_fold_change" FLOAT,
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
CREATE INDEX "isoformExpDiffData.fk_isoformExpDiffData_isoforms1" ON "isoformExpDiffData"("isoform_id");
CREATE INDEX "isoformExpDiffData.fk_isoformExpDiffData_samples1" ON "isoformExpDiffData"("sample_1");
CREATE INDEX "isoformExpDiffData.fk_isoformExpDiffData_samples2" ON "isoformExpDiffData"("sample_2");
DROP TABLE IF EXISTS "isoformFeatures";
CREATE TABLE "isoformFeatures"(
  "isoform_id" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_isoformFeatures_isoforms1"
    FOREIGN KEY("isoform_id")
    REFERENCES "isoforms"("isoform_id")
);
CREATE INDEX "isoformFeatures.fk_isoformFeatures_isoforms1" ON "isoformFeatures"("isoform_id");
COMMIT;
