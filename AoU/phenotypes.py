#!/usr/bin/env python3
import hail as hl
import pandas as pd
import numpy as np
import subprocess
import os, re
import importlib.util

from AoU.paths import *
from pathlib import Path

# in any pop, increased by >10x and with a count >30
VARIANT_BLACKLIST = ['chrM:12705:C,T', 'chrM:12684:G,A', 
                     'chrM:11467:A,G', 'chrM:13062:A,G', 'chrM:13095:T,C']

def get_final_annotated_variant_path(version):
    if version == 'v6':
        return os.path.join('gs://fc-secure-65229b17-5f6d-4315-9519-f53618eeee91/final_callset_220920/vcf/221012_filt_annotated/annotated_combined_processed_flat.tsv.bgz')
    if version == 'v7':
        return os.path.join(BUCKET, 'final_v7_merged_callset/all_v7_callset_20241030/vcf/annotated/annotated_combined_processed_flat.tsv.bgz')
    if version == 'v6andv7':
        raise NotImplementedError('It is VERY MUCH not recommended to use the combined callset for v6 + v7.')


def get_final_annotated_snv_flat_path(version, legacy=False):
    if version == 'v6':
        return os.path.join('gs://fc-secure-65229b17-5f6d-4315-9519-f53618eeee91/final_callset_220920/snv/221012_filt_annotated_all_heteroplasmic_snvs_case_only_v6_processed_flat.tsv.bgz')
    if version == 'v7':
        return os.path.join(BUCKET, 'final_v7_merged_callset/all_v7_callset_20241030/snv/annotated_all_heteroplasmic_snvs_case_only_v7_processed_flat.tsv.bgz')
    if version == 'v6andv7':
        if legacy:
            return os.path.join(BUCKET, 'full_250k_callset/snv/241128_filt_annotated_all_heteroplasmic_snvs_case_only_250k_processed_flat.tsv.bgz')
        else:
            return os.path.join(BUCKET, 'full_250k_callset/snv/250211_filt_annotated_1pct_all_heteroplasmic_snvs_case_only_250k_processed_flat.tsv.bgz')


def get_final_annotated_variant_mt_path(version, legacy=False):
    if version == 'v6':
        return os.path.join('gs://fc-secure-65229b17-5f6d-4315-9519-f53618eeee91/final_callset_220920/vcf/221012_filt_annotated/annotated_combined.mt')
    if version == 'v7':
        return os.path.join(BUCKET, 'final_v7_merged_callset/all_v7_callset_20241030/vcf/filt_annotated2/annotated_combined.mt')
    if version == 'v6andv7':
        if legacy:
            return os.path.join(BUCKET, 'full_250k_callset/vcf/241128_filt_annotated/annotated_combined.mt')
        else:
            return os.path.join(BUCKET, 'full_250k_callset/vcf/250211_filt_annotated_1pct/annotated_combined.mt')


def get_final_munged_case_only_hl_path(num_to_keep, type, version):
    base_path = os.path.join(BUCKET, f'munged_mtdna_callsets/{version}/munged/')
    if type == 'ht':
        return base_path + f'case_only_calls_low_hl_fullqc_N_{str(num_to_keep)}.ht'
    elif type == 'tsv':
        return base_path + f'case_only_calls_low_hl_fullqc_long_format_N_{str(num_to_keep)}.tsv'


def get_final_munged_snv_variant_path(version, HL):
    base_path = os.path.join(BUCKET, f'munged_mtdna_callsets/{version}/munged/')
    return base_path + f'heteroplasmic_snv_count_annotated_fullqc_{str(HL)}.mt'


def get_final_munged_snvcount_path(version, extension='ht'):
    base_path = os.path.join(BUCKET, f'munged_mtdna_callsets/{version}/munged/')
    return base_path + f'heteroplasmic_snv_indel_count_fullqc.{extension}'


def get_final_munged_snvcount_byclass_path(version, extension='ht'):
    base_path = os.path.join(BUCKET, f'munged_mtdna_callsets/{version}/munged/')
    return base_path + f'heteroplasmic_snv_count_by_class_fullqc.{extension}'


def get_final_munged_snvcount_age_accum_path(version, HL, extension='ht'):
    base_path = os.path.join(BUCKET, f'munged_mtdna_callsets/{version}/munged/')
    return base_path + f'heteroplasmic_snv_count_hl_age_accum_fullqc_{str(HL)}.{extension}'


def get_final_sample_stats_flat_path(version):
    if version == 'v6':
        return 'gs://fc-secure-65229b17-5f6d-4315-9519-f53618eeee91/final_callset_220920/221012_filtered_aou_tab_per_sample_stats.tsv'
    if version == 'v7':
        return os.path.join(BUCKET, 'final_v7_merged_callset/total_stats.merge.csv')
    if version == 'v6andv7':
        return os.path.join(BUCKET, 'full_250k_callset/241025_filtered_aou_total_stats.tsv')


def get_hap_ht_path(version, format, cutoff):
    base_path = os.path.join(BUCKET, f'munged_mtdna_callsets/{version}/munged/')
    if format == 'wide':
        return base_path + f'haplogroup_wide_N_{str(cutoff)}.ht'
    else:
        return base_path + f'haplogroup_tall.ht'


def get_hap_flat_path(version, format, cutoff):
    base_path = os.path.join(BUCKET, f'munged_mtdna_callsets/{version}/munged/')
    if format == 'wide':
        return base_path + f'haplogroup_wide_N_{str(cutoff)}.tsv.bgz'
    else:
        return base_path + f'haplogroup_tall.tsv.bgz'


def get_path_raw_positive_control():
    return os.path.join(PHENO_PATH, 'raw', '241016_positive_control_phenotypes.tsv')


def get_path_fancy_positive_control():
    return os.path.join(PHENO_PATH, 'raw/sanity_phenotypes', 'sanity_check_demographic_phenos.ht')


def get_path_fancy_positive_control_covars():
    return os.path.join(PHENO_PATH, 'raw/sanity_phenotypes', 'sanity_check_demographic_phenos_covars.ht')


def get_positive_control_phenotype_names():
    measure_list = {3036277: 'centimeter', 
                    3004501: 'milligram per deciliter', 
                    903115: 'millimeter mercury column',
                    903118: 'millimeter mercury column'}
    measures = '903133, 3004501, 903115, 903118'
    allowed_operators = ['=','No matching concept']
    return measure_list, measures, allowed_operators


def grab_positive_control_phenos(measures, measure_list, allowed_operators):
    data_sql = f"""
        SELECT 
            measurement.PERSON_ID,
            measurement.MEASUREMENT_DATETIME,
            m_ext.src_id as DATA_SOURCE,
            measurement.visit_occurrence_id as VISIT_ID,
            measurement.MEASUREMENT_CONCEPT_ID, 
            m_standard_concept.concept_name as STANDARD_CONCEPT_NAME, 
            m_standard_concept.vocabulary_id as STANDARD_VOCABULARY, 
            measurement.UNIT_CONCEPT_ID, 
            m_unit.concept_name as UNIT_CONCEPT_NAME, 
            measurement.MEASUREMENT_TYPE_CONCEPT_ID, 
            m_type.concept_name as MEASUREMENT_TYPE_CONCEPT_NAME, 
            measurement.MEASUREMENT_SOURCE_CONCEPT_ID, 
            m_source_concept.concept_name as SOURCE_CONCEPT_NAME, 
            m_source_concept.vocabulary_id as SOURCE_VOCABULARY,
            measurement.VALUE_AS_NUMBER,
            measurement.VALUE_AS_CONCEPT_ID, 
            m_value.concept_name as VALUE_AS_CONCEPT_NAME, 
            m_operator.concept_name as OPERATOR_CONCEPT_NAME
        FROM `{DATASET}.measurement` measurement 
        LEFT JOIN `{DATASET}.measurement_ext` m_ext 
            ON measurement.measurement_id = m_ext.measurement_id
        LEFT JOIN `{DATASET}.concept` m_standard_concept 
            ON measurement.measurement_concept_id = m_standard_concept.concept_id 
        LEFT JOIN `{DATASET}.concept` m_unit 
            ON measurement.unit_concept_id = m_unit.concept_id 
        LEFT JOIN `{DATASET}.concept` m_type 
            ON measurement.measurement_type_concept_id = m_type.concept_id 
        LEFT JOIN `{DATASET}.concept` m_source_concept 
            ON measurement.measurement_source_concept_id = m_source_concept.concept_id 
        LEFT JOIN `{DATASET}.concept` m_value 
            ON measurement.value_as_concept_id = m_value.concept_id 
        LEFT JOIN `{DATASET}.concept` m_operator 
            ON measurement.operator_concept_id = m_operator.concept_id 
        WHERE
            measurement_source_concept_id IN ({measures})"""
    
    wgs_data_measure = pd.read_gbq(data_sql, dialect="standard")
    wgs_data_measure_latest = wgs_data_measure.sort_values('MEASUREMENT_DATETIME'
                                                ).groupby(['MEASUREMENT_CONCEPT_ID', 'PERSON_ID']
                                                ).tail(1)
    wgs_data_measure_latest['correct_measure'] = wgs_data_measure_latest.MEASUREMENT_CONCEPT_ID.map(measure_list)
    wgs_data_measure_latest = wgs_data_measure_latest[wgs_data_measure_latest.correct_measure == wgs_data_measure_latest.UNIT_CONCEPT_NAME]
    wgs_data_measure_latest = wgs_data_measure_latest[wgs_data_measure_latest.OPERATOR_CONCEPT_NAME.isin(allowed_operators)]
    wgs_data_measure_latest['measurement_year'] = pd.DatetimeIndex(wgs_data_measure_latest.MEASUREMENT_DATETIME).year
    wgs_data_measure_latest['measurement_time'] = wgs_data_measure_latest.MEASUREMENT_DATETIME.astype('int') / 10**9

    return wgs_data_measure, wgs_data_measure_latest


def munge_phenotype_name(expr):
    """ Applies basic transformations to make phenotype names legible
    """
    return expr.replace('[\\[\\]\\(\\)/\s,.\\\\:\'\"]{1,}', '_').lower()


def import_and_process_pos_control_phenos(path_latest_csv):
    ht_latest = hl.import_table(path_latest_csv, min_partitions=100, impute=True, 
                                types={'VALUE_AS_NUMBER': hl.tfloat64}, missing=['NA',''])
    ht_latest = ht_latest.annotate(phenotype = munge_phenotype_name(ht_latest.STANDARD_CONCEPT_NAME))
    ht_latest = ht_latest.annotate(s = hl.str(ht_latest.PERSON_ID)).drop('PERSON_ID').key_by('s')
    return ht_latest


def get_raw_positive_control_phenotypes():
    pheno_path = get_path_raw_positive_control()
    if not hl.hadoop_exists(pheno_path):
        measure_list, measures, allowed_operators = get_positive_control_phenotype_names()
        _, wgs_data_measure_latest = grab_positive_control_phenos(measures=measures, measure_list=measure_list, allowed_operators=allowed_operators)
        
        wgs_data_measure_latest['munged_concept_name'] = wgs_data_measure_latest.STANDARD_CONCEPT_NAME.str.replace('[\\[\\]\\(\\)/\s,.\\\\:\'\"]{1,}', '_').str.lower()
        backbone = wgs_data_measure_latest[['PERSON_ID']].drop_duplicates().rename({'PERSON_ID': 's'}, axis=1)
        phenotypes = set(wgs_data_measure_latest['munged_concept_name'])
        pheno_map = {'body height': 'height',
                    'computed diastolic blood pressure, mean of 2nd and 3rd measures': 'diastolic_bp',
                    'computed systolic blood pressure, mean of 2nd and 3rd measures': 'systolic_bp',
                    'glucose [mass/volume] in serum or plasma': 'glucose'}
        
        for pheno in phenotypes:
            this_target_name = pheno_map[pheno]
            df_sub = wgs_data_measure_latest[wgs_data_measure_latest.munged_concept_name == pheno][['PERSON_ID', 'VALUE_AS_NUMBER']].rename({'PERSON_ID': 's', 'VALUE_AS_NUMBER': this_target_name}, axis=1)
            backbone = backbone.merge(df_sub, how='left', on=['s'])

        local_pheno_path = os.path.abspath('./241016_positive_control_phenotypes.tsv')
        backbone.to_csv(local_pheno_path, sep='\t', index=False, na_rep='NA')
        _ = subprocess.run(['gsutil','-u',os.getenv('GOOGLE_PROJECT'),'cp',local_pheno_path, pheno_path])
    
    return backbone


def get_fancy_positive_control_phenotypes(sample_covariates, overwrite=False):

    measure_list, measures, allowed_operators = get_positive_control_phenotype_names()
    sanity_ht_path = get_path_fancy_positive_control()
    covar_ht_path = get_path_fancy_positive_control_covars()
    if hl.hadoop_exists(f'{sanity_ht_path}/_SUCCESS') and not overwrite:
        ht_pheno_out = hl.read_table(sanity_ht_path)
        ht_covar_out = hl.read_table(covar_ht_path)
        phenotypes = [x for x in ht_pheno_out.row if x not in ht_pheno_out.key]
    else:
        temp_phenos = f'{PHENO_PATH}/raw/sanity_phenotypes/demog_obtained_traits.tsv'
        temp_phenos_latest = f'{PHENO_PATH}/raw/sanity_phenotypes/demog_obtained_traits_latest.tsv'
        
        if not hl.hadoop_exists(temp_phenos) or overwrite:
            wgs_data_measure, wgs_data_measure_latest = grab_positive_control_phenos(measure_list=measure_list, measures=measures, allowed_operators=allowed_operators)
            wgs_data_measure.to_csv(temp_phenos, index=False, sep='\t')
            wgs_data_measure_latest.to_csv(temp_phenos_latest, index=False, sep='\t')

        ht_latest = import_and_process_pos_control_phenos(temp_phenos_latest)
        ht_latest = ht_latest.annotate(year_of_birth = sample_covariates[ht_latest.s].year_of_birth,
                                       isFemale = sample_covariates[ht_latest.s].isFemale)
        ht_latest = ht_latest.annotate(approx_age = ht_latest.measurement_year - ht_latest.year_of_birth)
        ht_latest = ht_latest.filter(ht_latest.approx_age > 0)
        ht_latest = ht_latest.checkpoint(f'{TEMP_PATH}/sanity_check_alldata_checkpoint.ht', overwrite=True)
        ht_pheno_out = ht_latest.select().distinct()
        ht_covar_out = ht_latest.select().distinct()
        phenotypes = list(ht_latest.aggregate(hl.agg.collect_as_set(ht_latest.phenotype)))

        for pheno in phenotypes:
            htf = ht_latest.filter(ht_latest.phenotype == pheno)
            ht_pheno_out = ht_pheno_out.annotate(**{pheno: htf[ht_pheno_out.s].VALUE_AS_NUMBER})

            struct_item = hl.Struct(**{'measurement_time': htf[ht_covar_out.s].measurement_time,
                                       'collection_site': htf[ht_covar_out.s].DATA_SOURCE,
                                       'visit_id': htf[ht_covar_out.s].VISIT_ID,
                                       'approx_age': htf[ht_covar_out.s].approx_age,
                                       'age2': htf[ht_covar_out.s].approx_age**2,
                                       'age_isFemale':htf[ht_covar_out.s].isFemale * htf[ht_covar_out.s].approx_age,
                                       'age2_isFemale': (htf[ht_covar_out.s].approx_age**2)*htf[ht_covar_out.s].isFemale})
            ht_covar_out = ht_covar_out.annotate(**{pheno: struct_item})
            ht_covar_out = ht_covar_out.checkpoint(f'{TEMP_PATH}/sanity_check_covars_after_{pheno}_checkpoint.ht', overwrite=True)

        ht_pheno_out = ht_pheno_out.checkpoint(sanity_ht_path, overwrite=overwrite)
        ht_covar_out = ht_covar_out.checkpoint(covar_ht_path, overwrite=overwrite)

    return ht_pheno_out, ht_covar_out, phenotypes


def get_case_only_mtdna_callset(version, num_to_keep=500, overwrite=False):
    """ Here we implement the munging pipeline in Python and Pandas.
    """
    if hl.hadoop_exists(f"{get_final_munged_case_only_hl_path(num_to_keep, 'ht', version)}/_SUCCESS") and not overwrite:
        ht_wide_hl_all = hl.read_table(get_final_munged_case_only_hl_path(num_to_keep, 'ht', version))
    
    else:
        mt = hl.read_matrix_table(get_final_annotated_variant_mt_path(version))
        ht_N_filt = hl.import_table(f'{BUCKET}reference_files/heteroplasmic_variants_ukb_aou_{str(num_to_keep)}.tsv',
                                    types={'locus': hl.tlocus(reference_genome='GRCh38'),
                                           'alleles': hl.tarray('str')}, impute=True)
        ht_N_filt = ht_N_filt.key_by('locus','alleles')
        mt_variants_hl = mt.semi_join_rows(ht_N_filt).select_rows().select_cols().select_globals()
        mt_variants_hl = mt_variants_hl.annotate_entries(HL = hl.case().when((mt_variants_hl.HL >= 0.05) & (mt_variants_hl.HL < 0.95) & (hl.is_defined(mt_variants_hl.HL)), mt_variants_hl.HL
                                                                    ).or_missing())
        ht_variants_hl = mt_variants_hl.select_entries('HL').entries()
        ht_variants_hl = ht_variants_hl.annotate(variant = ht_variants_hl.locus.contig + '_' + hl.str(ht_variants_hl.locus.position) + '_' + hl.str('_').join(ht_variants_hl.alleles))
        ht_variants_hl = ht_variants_hl.key_by('s').select('variant', 'HL').checkpoint(f'{TEMP_PATH}/variants_entries_for_low_HL.ht', overwrite=True)
        variants = ht_variants_hl.aggregate(hl.agg.collect_as_set(ht_variants_hl.variant))
        variants = list(np.array(list(variants))[np.argsort([int(x.split('_')[1]) for x in variants])])
        iter_size = 10
        hl_list = []
        for idx, start in enumerate(range(0, len(variants), iter_size)):
            end = start + iter_size if start + iter_size < len(variants) else len(variants)
            these_vars = variants[start:end]
            print('Making the following variants wide:')
            print(these_vars)
            ht_wide_hl = mt_variants_hl.cols().select()
            ht_wide_hl = ht_wide_hl.annotate(**{x: ht_variants_hl.filter(ht_variants_hl.variant == x)[ht_wide_hl.key].HL for x in these_vars})
            ht_wide_hl = ht_wide_hl.checkpoint(f'{TEMP_PATH}/variants_hl_temp_{str(idx)}.ht', overwrite=True)
            hl_list.append(ht_wide_hl)
        ht_wide_hl_all = hl_list[0]
        for this_ht in hl_list[1:len(hl_list)]:
            ht_wide_hl_all = ht_wide_hl_all.join(this_ht, how='outer')
        ht_wide_hl_all = ht_wide_hl_all.checkpoint(get_final_munged_case_only_hl_path(num_to_keep, 'ht', version), overwrite=True)
        ht_wide_hl_all.export(get_final_munged_case_only_hl_path(num_to_keep, 'tsv', version))
    
    return ht_wide_hl_all


def get_snv_indel_count_phenotype(version, overwrite=False):
    if hl.hadoop_exists(f"{get_final_munged_snvcount_path(version)}/_SUCCESS") and not overwrite:
        ht_snv_count = hl.read_table(get_final_munged_snvcount_path(version))
    
    else:
        mt = hl.read_matrix_table(get_final_annotated_variant_mt_path(version))
        mt_snv = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1])).select_cols()
        mt_snv = mt_snv.annotate_cols(snv_count = hl.agg.count_where(hl.is_defined(mt_snv.HL) & (mt_snv.HL < 0.95) & (mt_snv.HL >= 0.05)))
        ht_snv_count = mt_snv.cols()

        mt_indel = mt.filter_rows(hl.is_indel(mt.alleles[0], mt.alleles[1])).select_cols()
        mt_indel = mt_indel.annotate_cols(indel_count = hl.agg.count_where(hl.is_defined(mt_indel.HL) & (mt_indel.HL < 0.95) & (mt_indel.HL >= 0.05)))
        ht_count = mt_indel.cols()
        ht_count = ht_count.annotate(snv_count = ht_snv_count[ht_count.key].snv_count)

        common_indel_sites = [302, 310, 513, 567, 955, 5894, 8270, 16182, 16183, 16189, 16192]
        mt_indel_uncommon = mt.filter_rows(hl.is_indel(mt.alleles[0], mt.alleles[1])).select_cols()
        mt_indel_uncommon = mt_indel_uncommon.filter_rows(~hl.literal(common_indel_sites).contains(mt_indel_uncommon.locus.position))
        mt_indel_uncommon = mt_indel_uncommon.annotate_cols(indel_count = hl.agg.count_where(hl.is_defined(mt_indel_uncommon.HL) & (mt_indel_uncommon.HL < 0.95) & (mt_indel_uncommon.HL >= 0.05)))
        ht_count = ht_count.annotate(indel_count_uncommon = mt_indel_uncommon.cols()[ht_count.key].indel_count)

        mt_indel_common = mt.filter_rows(hl.is_indel(mt.alleles[0], mt.alleles[1])).select_cols()
        mt_indel_common = mt_indel_common.filter_rows(hl.literal(common_indel_sites).contains(mt_indel_common.locus.position))
        mt_indel_common = mt_indel_common.annotate_cols(indel_count = hl.agg.count_where(hl.is_defined(mt_indel_common.HL) & (mt_indel_common.HL < 0.95) & (mt_indel_common.HL >= 0.05)))
        ht_count = ht_count.annotate(indel_count_common = mt_indel_common.cols()[ht_count.key].indel_count)

        ht_snv_count = ht_count.select_globals().checkpoint(get_final_munged_snvcount_path(version), overwrite=True)
        ht_snv_count.export(get_final_munged_snvcount_path(version, extension='tsv'))
    
    return ht_snv_count


def get_snv_count_by_class(version, overwrite=False):
    if hl.hadoop_exists(f"{get_final_munged_snvcount_byclass_path(version)}/_SUCCESS") and not overwrite:
        ht_wide_snv_count = hl.read_table(get_final_munged_snvcount_byclass_path(version))
    
    else:
        mt = hl.read_matrix_table(get_final_annotated_variant_mt_path(version))
        mt_snv_class = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1])).select_cols()
        mt_snv_class = mt_snv_class.annotate_rows(allele_string = hl.str('_').join(mt_snv_class.alleles))
        mt_snv_count_by_class = mt_snv_class.group_rows_by(mt_snv_class.allele_string).aggregate(snv_count = hl.agg.count_where(hl.is_defined(mt_snv_class.HL) & (mt_snv_class.HL < 0.95) & (mt_snv_class.HL >= 0.05))).select_globals()

        ht_snv_count_by_class = mt_snv_count_by_class.entries()
        ht_snv_count_by_class = ht_snv_count_by_class.checkpoint(TEMP_PATH + '/snv_count_by_class_long_intermediate.ht', overwrite=True)
        all_classes = list(ht_snv_count_by_class.aggregate(hl.agg.collect_as_set(ht_snv_count_by_class.allele_string)))
        this_ht_snv_count_by_class = ht_snv_count_by_class.key_by('s')
        ht_wide_snv_count = mt_snv_count_by_class.cols().select()
        ht_wide_snv_count = ht_wide_snv_count.annotate(**{f'snv_count_{x}': this_ht_snv_count_by_class.filter(this_ht_snv_count_by_class.allele_string == x)[ht_wide_snv_count.key].snv_count for x in all_classes})
        ht_wide_snv_count = ht_wide_snv_count.checkpoint(get_final_munged_snvcount_byclass_path(version), overwrite=True)
        ht_wide_snv_count.export(get_final_munged_snvcount_byclass_path(version, 'tsv'))
    
    return ht_wide_snv_count


def get_annotated_snvs(version, HL, overwrite=False):
    if hl.hadoop_exists(f'{get_final_munged_snv_variant_path(version, HL)}/_SUCCESS') and not overwrite:
        mt_snv_class = hl.read_matrix_table(get_final_munged_snv_variant_path(version, HL))
    
    else:
        mt = hl.read_matrix_table(get_final_annotated_variant_mt_path('v6andv7'))
        bp_compl = hl.literal({'C': 'G', 'G': 'C', 'T': 'A', 'A': 'T'})

        # add variant class and location annotations
        mt_snv_class = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1])).select_cols()
        mt_snv_class = mt_snv_class.annotate_rows(allele_compl = [bp_compl[mt_snv_class.alleles[0]], bp_compl[mt_snv_class.alleles[1]]])
        mt_snv_class = mt_snv_class.annotate_rows(location = hl.if_else((mt_snv_class.locus.position > 16172) | (mt_snv_class.locus.position < 210), 'Ori (16172-210)', 'Other region'))
        mt_snv_class = mt_snv_class.annotate_rows(variant_class = hl.if_else(hl.literal(['C','T']).contains(mt_snv_class.alleles[0]), 
                                                                            mt_snv_class.alleles[0] + '>' + mt_snv_class.alleles[1],
                                                                            mt_snv_class.allele_compl[0] + '>' + mt_snv_class.allele_compl[1]),
                                                  strand = hl.if_else(hl.literal(['C','T']).contains(mt_snv_class.alleles[0]), 'light', 'heavy'))
        mt_snv_class = mt_snv_class.drop('allele_compl')
        mt_snv_class = mt_snv_class.annotate_rows(variant = mt_snv_class.locus.contig + ':' + \
                                                                hl.str(mt_snv_class.locus.position) + ':' + \
                                                                hl.str(',').join(mt_snv_class.alleles))

        # add poisson-based cutoffs
        ht_stats = hl.import_table(get_final_sample_stats_flat_path(version=version), impute=True).select('s','nuc_mean_coverage')
        ht_stats = ht_stats.annotate(s = hl.str(ht_stats.s))
        ht_stats = ht_stats.transmute(pois_cutoff = hl.qpois(0.95, ht_stats.nuc_mean_coverage),
                                      pois_cutoff975 = hl.qpois(0.975, ht_stats.nuc_mean_coverage),
                                      pois_cutoff99 = hl.qpois(0.99, ht_stats.nuc_mean_coverage)).key_by('s')
        mt_snv_class = mt_snv_class.annotate_cols(pois_cutoff = ht_stats[mt_snv_class.col_key].pois_cutoff,
                                                  pois_cutoff975 = ht_stats[mt_snv_class.col_key].pois_cutoff975,
                                                  pois_cutoff99 = ht_stats[mt_snv_class.col_key].pois_cutoff99)
        
        # save
        mt_snv_class = mt_snv_class.checkpoint(get_final_munged_snv_variant_path(version, HL), overwrite=True)
    
    return mt_snv_class


def get_age_accumulating_snv_count(version, HL, overwrite=False):
    if hl.hadoop_exists(f"{get_final_munged_snvcount_age_accum_path(version, HL)}/_SUCCESS") and not overwrite:
        ht_snv_age_accum_count_final = hl.read_table(get_final_munged_snvcount_age_accum_path(version, HL))
    
    else:
        mt_snv_class = get_annotated_snvs(version, HL, overwrite)
        mt_snv_class_f = mt_snv_class.filter_rows((mt_snv_class.variant_class == 'T>C') | ((mt_snv_class.variant_class == 'C>T') & (mt_snv_class.strand == 'heavy')))
        mt_snv_class_f = mt_snv_class_f.filter_rows(mt_snv_class_f.location != 'Ori (16172-210)')
        ht_snv_age_accum = mt_snv_class_f.entries()

        # sanity checks
        # ht_this = ht_snv_age_accum.filter(ht_snv_age_accum.AD[1] > ht_snv_age_accum.pois_cutoff)
        # ht_this = ht_this.filter(hl.is_defined(ht_this.HL) & (ht_this.HL < 0.95) & (ht_this.HL >= 0.01))
        # ht_this.group_by(ht_this.variant).aggregate(ct = hl.agg.count()).order_by(hl.desc('ct')).show(200)
        # | "chrM:12684:G,A" |  7099 | <- on the variant blacklist!
        # | "chrM:16093:T,C" |  3120 |
        # | "chrM:204:T,C"   |  2059 |
        # | "chrM:16129:G,A" |  1964 |
        # | "chrM:146:T,C"   |  1726 |
        # | "chrM:152:T,C"   |  1493 |
        # | "chrM:16311:T,C" |  1128 |
        # | "chrM:195:T,C"   |   757 |
        # | "chrM:13062:A,G" |   734 | <- blacklisted
        # | "chrM:16519:T,C" |   610 |
        # | "chrM:2623:A,G"  |   599 |
        # | "chrM:16189:T,C" |   512 |
        # | "chrM:16362:T,C" |   510 |


        def _generate_all_traits(ht, HL):
            tf_expression = hl.is_defined(ht.HL) & (ht.HL < 0.95)
            ht_count = ht.group_by(ht.s
                        ).aggregate(snv_count = hl.agg.count_where(tf_expression & (ht.HL >= HL)),
                                    snv_invHL = hl.agg.filter(tf_expression & (ht.HL >= HL), hl.agg.sum(1 / ht.HL)),
                                    snv_sumHL = hl.agg.filter(tf_expression & (ht.HL >= HL), hl.agg.sum(ht.HL)),
                                    snv_1minHL = hl.agg.filter(tf_expression & (ht.HL >= HL), hl.agg.sum(1 - ht.HL)),
                                    snv_meanHL = hl.agg.filter(tf_expression & (ht.HL >= HL), hl.agg.mean(ht.HL)),
                                    snv_ADpois95_count = hl.agg.count_where(tf_expression & (ht.HL >= 0.01) & (ht.AD[1] > ht.pois_cutoff)),
                                    snv_ADpois95_invHL = hl.agg.filter(tf_expression & (ht.HL >= 0.01) & (ht.AD[1] > ht.pois_cutoff), hl.agg.sum(1 / ht.HL)),
                                    snv_ADpois95_1minHL = hl.agg.filter(tf_expression & (ht.HL >= 0.01) & (ht.AD[1] > ht.pois_cutoff), hl.agg.sum(1 - ht.HL)),
                                    snv_ADpois975_1minHL = hl.agg.filter(tf_expression & (ht.HL >= 0.01) & (ht.AD[1] > ht.pois_cutoff975), hl.agg.sum(1 - ht.HL)),
                                    snv_ADpois975_count = hl.agg.count_where(tf_expression & (ht.HL >= 0.01) & (ht.AD[1] > ht.pois_cutoff975)),
                                    snv_ADpois99_1minHL = hl.agg.filter(tf_expression & (ht.HL >= 0.01) & (ht.AD[1] > ht.pois_cutoff99), hl.agg.sum(1 - ht.HL)),
                                    snv_ADpois99_count = hl.agg.count_where(tf_expression & (ht.HL >= 0.01) & (ht.AD[1] > ht.pois_cutoff99))).select_globals()
            ht_count = ht_count.annotate(snv_meanHL = hl.if_else(hl.is_nan(ht_count.snv_meanHL), 0, ht_count.snv_meanHL))
            return ht_count


        print('Generating phenotypes with and without blacklisted sites...')
        print('PLEASE NOTE: a ton of phenotypes will be generated. Please subset to the desired traits for analysis.')
        ht_snv_age_accum_count = _generate_all_traits(ht_snv_age_accum, HL=HL)
        ht_snv_age_accum_count = ht_snv_age_accum_count.rename({x: f'{x}_nosnvrm' for x in ht_snv_age_accum_count.row if x not in ht_snv_age_accum_count.key})
        ht_snv_age_accum_count = ht_snv_age_accum_count.checkpoint(os.path.join(TEMP_PATH, 'ht_snv_age.ht'), overwrite=True)

        ht_snv_age_accum_f = ht_snv_age_accum.filter(~hl.literal(VARIANT_BLACKLIST).contains(ht_snv_age_accum.variant))
        ht_snv_age_accum_count_f = _generate_all_traits(ht_snv_age_accum_f, HL=HL)
        ht_snv_age_accum_count_f = ht_snv_age_accum_count_f.checkpoint(os.path.join(TEMP_PATH, 'ht_snv_age_filt.ht'), overwrite=True)

        print('Generating phenotypes without CHIP individuals...')
        chip_data = load_chip_data()
        chip_data = chip_data.filter(chip_data.hasCH != 0)
        ht_snv_age_accum_nochip = mt_snv_class_f.anti_join_cols(chip_data).entries()
        ht_snv_age_accum_f_nochip = ht_snv_age_accum_nochip.filter(~hl.literal(VARIANT_BLACKLIST).contains(ht_snv_age_accum_nochip.variant))
        ht_snv_age_accum_count_f_nochip = _generate_all_traits(ht_snv_age_accum_f_nochip, HL=HL)
        ht_snv_age_accum_count_f_nochip = ht_snv_age_accum_count_f_nochip.rename({x: f'{x}_nochip' for x in ht_snv_age_accum_count_f_nochip.row if x not in ht_snv_age_accum_count_f_nochip.key})
        ht_snv_age_accum_count_f_nochip = ht_snv_age_accum_count_f_nochip.checkpoint(os.path.join(TEMP_PATH, 'ht_snv_age_filt_nochip.ht'), overwrite=True)

        print('Combining all traits...')
        ht_snv_age_accum_count_final = ht_snv_age_accum_count_f.annotate(**ht_snv_age_accum_count[ht_snv_age_accum_count_f.key],
                                                                         **ht_snv_age_accum_count_f_nochip[ht_snv_age_accum_count_f.key])
        ht_snv_age_accum_count_final = ht_snv_age_accum_count_final.checkpoint(get_final_munged_snvcount_age_accum_path(version, HL), overwrite=True)
        ht_snv_age_accum_count_final.export(get_final_munged_snvcount_age_accum_path(version, HL, 'tsv'))
    
    return ht_snv_age_accum_count_final


def extract_single_mtdna_variant(position, ref, alt, version):
    mt = hl.read_matrix_table(get_final_annotated_variant_mt_path(version))
    mtf = mt.filter_rows((mt.locus.position == position) & (mt.alleles == [ref, alt]))

    mtf = mtf.select_globals().select_cols()
    mtf = mtf.key_rows_by(variant = mtf.locus.contig + ':' + hl.str(mtf.locus.position) + ':' + hl.str(':').join(mtf.alleles)).select_rows('filters')
    ht = mtf.entries()
    ht = ht.annotate(AD_ref = ht.AD[0], AD_alt = ht.AD[1], FT = ht.FT.union(ht.filters)).drop('AD','filters')
    ht = ht.annotate(OriginalSelfRefAlleles = hl.str(',').join(ht.OriginalSelfRefAlleles))
    if 'F1R2' in ht.row and 'F2R1' in ht.row:
        ht = ht.annotate(F2R1_ref = ht.F2R1[0], F2R1_alt = ht.F2R1[1], F1R2_ref = ht.F1R2[0], F1R2_alt = ht.F1R2[1]).drop('F2R1','F1R2')
    if 'RPA' in ht.row:
        ht = ht.annotate(RPA_ref = ht.RPA[0], RPA_alt = ht.RPA[1]).drop('RPA')
    if 'AS_SB_TABLE' in ht.row:
        ht = ht.annotate(AS_SB_SPLIT = ht.AS_SB_TABLE.split('\\|'))
        ht = ht.annotate(FWD_SB_ref = hl.if_else(hl.is_defined(ht.AS_SB_TABLE), hl.int32(ht.AS_SB_SPLIT[0].split(',')[0]), hl.missing(hl.tint32)),
                            FWD_SB_alt = hl.if_else(hl.is_defined(ht.AS_SB_TABLE), hl.int32(ht.AS_SB_SPLIT[1].split(',')[0]), hl.missing(hl.tint32)),
                            REV_SB_ref = hl.if_else(hl.is_defined(ht.AS_SB_TABLE), hl.int32(ht.AS_SB_SPLIT[0].split(',')[1]), hl.missing(hl.tint32)),
                            REV_SB_alt = hl.if_else(hl.is_defined(ht.AS_SB_TABLE), hl.int32(ht.AS_SB_SPLIT[1].split(',')[1]), hl.missing(hl.tint32)))
        ht = ht.drop('AS_SB_TABLE', 'AS_SB_SPLIT')
    ht = ht.annotate(FT = hl.literal(',').join(ht.FT.difference({'PASS'})), 
                    FT_LIFT = hl.literal(',').join(ht.FT_LIFT.difference({'PASS'}))).drop('GT')

    return ht


def load_chip_data():
    CHIP_path = os.path.join(BUCKET, 'CHIP_240902/traits/')
    CHIP_ht = os.path.join(CHIP_path, 'munged/aou_v7_chip_calls.ht')

    if hl.hadoop_exists(os.path.join(CHIP_ht, '_SUCCESS')):
        ht = hl.read_table(CHIP_ht)
    else:
        aou_chip_all_var = hl.import_table(os.path.join(CHIP_path, 'aou_98k_150k.chip_var.dp20ad3ad5.allVar.26Aug2023.csv.gz'), 
                                           types={'SampleID': hl.tstr}, delimiter=',', impute=True, force=True, quote='"')
        aou_chip_pheno = hl.import_table(os.path.join(CHIP_path, 'phenoCH_AoU250k.fid0_iid.qcd_myeloidCA_rel_NA.26Aug2023.tsv.gz'), 
                                         types={'IID': hl.tstr}, impute=True, force=True)
        print('AoU CHIP data imported.')
        
        aou_top_chip_table = aou_chip_all_var.group_by(aou_chip_all_var.SampleID).aggregate(VAF = hl.agg.max(aou_chip_all_var.VAF), batch = hl.agg.take(aou_chip_all_var.Batch_CH, 1))
        aou_jt_chip_table = aou_top_chip_table.key_by('SampleID', 'VAF').join(aou_chip_all_var.select('SampleID','varID','Gene','Nonsyn','DP','AD','VAF').key_by('SampleID','VAF'), how='left')
        aou_CHIP_calls = aou_chip_pheno.select(s = aou_chip_pheno.IID, age = aou_chip_pheno.Age_biosample_collection, hasCH = aou_chip_pheno.hasCH, hasCHvaf10 = aou_chip_pheno.hasCHvaf10)
        aou_final_chip_calls = aou_CHIP_calls.key_by('s').join(aou_jt_chip_table.key_by('SampleID'), how='left')
        aou_final_chip_calls = aou_final_chip_calls.drop('age')
        aou_final_chip_calls = aou_final_chip_calls.annotate(dataset = 'AoU')

        aou_dupes = aou_final_chip_calls.group_by(aou_final_chip_calls.s).aggregate(N = hl.agg.count())
        aou_dupes = aou_dupes.filter(aou_dupes.N > 1)
        aou_final_chip_calls_singular = aou_final_chip_calls.anti_join(aou_dupes)
        aou_final_chip_calls_duplicates = aou_final_chip_calls.semi_join(aou_dupes).add_index()
        aou_final_chip_calls_idx = aou_final_chip_calls_duplicates.group_by(aou_final_chip_calls_duplicates.s, aou_final_chip_calls_duplicates.VAF).aggregate(idx = hl.agg.collect_as_set(aou_final_chip_calls_duplicates.idx))
        aou_final_chip_calls_idx = aou_final_chip_calls_idx.annotate(idx = hl.sorted(aou_final_chip_calls_idx.idx))
        aou_final_chip_calls_idx = aou_final_chip_calls_idx.annotate(idx_flat = hl.find(lambda x: True, aou_final_chip_calls_idx.idx)).key_by('idx_flat')
        aou_final_chip_calls_dedupe = aou_final_chip_calls_duplicates.key_by('idx').semi_join(aou_final_chip_calls_idx).key_by('s').drop('idx')
        aou_final_chip_calls_proc = hl.Table.union(aou_final_chip_calls_singular, aou_final_chip_calls_dedupe)
        ht = aou_final_chip_calls_proc.checkpoint(CHIP_ht)
        ht.export(os.path.splitext(CHIP_ht)[0] + '.tsv.bgz')
    
    return ht


def export_all_snvs_to_flat_file(module_path, version, legacy=False):

    # Load in the processor function. Expected to be obtained from mtSwirl add_annotations.py.
    file_path = Path(module_path)
    module_name = file_path.stem
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    process_mt_for_flat_file_analysis = getattr(module, 'process_mt_for_flat_file_analysis')
    
    # load the callset
    mt = hl.read_matrix_table(get_final_annotated_variant_mt_path(version, legacy=legacy))
    mt_snv = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))
    mt_snv = mt_snv.annotate_rows(n_heteroplasmic = hl.agg.count_where(hl.is_defined(mt_snv.HL) & (mt_snv.HL < 0.95) & (mt_snv.HL >= 0.05)))
    mt_snv = mt_snv.filter_rows(mt_snv.n_heteroplasmic > 0)
    mt_snv.count() # (12973, 237500)

    mt_snv_export = mt_snv.drop('n_heteroplasmic')
    ht_export, _, _ = process_mt_for_flat_file_analysis(mt_snv_export, skip_vep=True, allow_gt_fail=False)
    ht_export = ht_export.filter(hl.is_defined(ht_export.HL) & (ht_export.HL < 0.95) & (ht_export.HL >= 0.01)) # 625708
    ht_export.export(get_final_annotated_snv_flat_path(version, legacy))