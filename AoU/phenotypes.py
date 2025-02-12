#!/usr/bin/env python3
import hail as hl
import pandas as pd
import numpy as np
import subprocess
import os, re

from AoU.paths import *
VARIANT_BLACKLIST = ['chrM:12705:C,T', 'chrM:12684:G,A']

def get_final_annotated_variant_path(version):
    if version == 'v6':
        return os.path.join('gs://fc-secure-65229b17-5f6d-4315-9519-f53618eeee91/final_callset_220920/vcf/221012_filt_annotated/annotated_combined_processed_flat.tsv.bgz')
    if version == 'v7':
        return os.path.join(BUCKET, 'final_v7_merged_callset/all_v7_callset_20241030/vcf/annotated/annotated_combined_processed_flat.tsv.bgz')
    if version == 'v6andv7':
        raise NotImplementedError('It is VERY MUCH not recommended to use the combined callset for v6 + v7.')


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


def get_age_accumulating_snv_count(version, HL, overwrite=False):
    if hl.hadoop_exists(f"{get_final_munged_snvcount_age_accum_path(version, HL)}/_SUCCESS") and not overwrite:
        ht_snv_age_accum_count = hl.read_table(get_final_munged_snvcount_age_accum_path(version, HL))
    
    else:
        mt = hl.read_matrix_table(get_final_annotated_variant_mt_path('v6andv7'))
        mt_snv_class = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1])).select_cols()

        ht_stats = hl.import_table(get_final_sample_stats_flat_path(version=version), impute=True).select('s','nuc_mean_coverage')
        ht_stats = ht_stats.annotate(s = hl.str(ht_stats.s))
        ht_stats = ht_stats.transmute(pois_cutoff = hl.qpois(0.95, ht_stats.nuc_mean_coverage)).key_by('s')

        bp_compl = hl.literal({'C': 'G', 'G': 'C', 'T': 'A', 'A': 'T'})
        mt_snv_class = mt_snv_class.annotate_rows(allele_compl = [bp_compl[mt_snv_class.alleles[0]], bp_compl[mt_snv_class.alleles[1]]])
        mt_snv_class = mt_snv_class.annotate_rows(location = hl.if_else((mt_snv_class.locus.position > 16172) | (mt_snv_class.locus.position < 210), 'Ori (16172-200)', 'Other region'))
        mt_snv_class = mt_snv_class.annotate_rows(variant_class = hl.if_else(hl.literal(['C','T']).contains(mt_snv_class.alleles[0]), 
                                                                            mt_snv_class.alleles[0] + '>' + mt_snv_class.alleles[1],
                                                                            mt_snv_class.allele_compl[0] + '>' + mt_snv_class.allele_compl[1]),
                                                  strand = hl.if_else(hl.literal(['C','T']).contains(mt_snv_class.alleles[0]), 'light', 'heavy'))
        mt_snv_class_f = mt_snv_class.filter_rows((mt_snv_class.variant_class == 'T>C') | ((mt_snv_class.variant_class == 'C>T') & (mt_snv_class.strand == 'heavy')))
        mt_snv_class_f = mt_snv_class_f.drop('allele_compl')
        mt_snv_class_f = mt_snv_class_f.annotate_cols(pois_cutoff = ht_stats[mt_snv_class_f.col_key].pois_cutoff)

        ht_snv_age_accum = mt_snv_class_f.entries()
        ht_snv_age_accum = ht_snv_age_accum.annotate(variant = ht_snv_age_accum.locus.contig + ':' + hl.str(ht_snv_age_accum.locus.position) + ':' + hl.str(',').join(ht_snv_age_accum.alleles))

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
        # | "chrM:13062:A,G" |   734 |
        # | "chrM:16519:T,C" |   610 |
        # | "chrM:2623:A,G"  |   599 |
        # | "chrM:16189:T,C" |   512 |
        # | "chrM:16362:T,C" |   510 |

        ht_snv_age_accum = ht_snv_age_accum.filter(~hl.literal(VARIANT_BLACKLIST).contains(ht_snv_age_accum.variant))
        ht_snv_age_accum_count = ht_snv_age_accum.group_by(ht_snv_age_accum.s
                                                ).aggregate(snv_count = hl.agg.count_where(hl.is_defined(ht_snv_age_accum.HL) & (ht_snv_age_accum.HL < 0.95) & (ht_snv_age_accum.HL >= HL)),
                                                            snv_invHL = hl.agg.filter(hl.is_defined(ht_snv_age_accum.HL) & (ht_snv_age_accum.HL < 0.95) & (ht_snv_age_accum.HL >= HL), hl.agg.sum(1 / ht_snv_age_accum.HL)),
                                                            snv_sumHL = hl.agg.filter(hl.is_defined(ht_snv_age_accum.HL) & (ht_snv_age_accum.HL < 0.95) & (ht_snv_age_accum.HL >= HL), hl.agg.sum(ht_snv_age_accum.HL)),
                                                            snv_1minHL = hl.agg.filter(hl.is_defined(ht_snv_age_accum.HL) & (ht_snv_age_accum.HL < 0.95) & (ht_snv_age_accum.HL >= HL), hl.agg.sum(1 - ht_snv_age_accum.HL)),
                                                            snv_meanHL = hl.agg.filter(hl.is_defined(ht_snv_age_accum.HL) & (ht_snv_age_accum.HL < 0.95) & (ht_snv_age_accum.HL >= HL), hl.agg.mean(ht_snv_age_accum.HL)),
                                                            snv_ADpois95_count = hl.agg.count_where(hl.is_defined(ht_snv_age_accum.HL) & (ht_snv_age_accum.HL < 0.95) & (ht_snv_age_accum.HL >= 0.01) & (ht_snv_age_accum.AD[1] > ht_snv_age_accum.pois_cutoff)),
                                                            snv_ADpois95_invHL = hl.agg.filter(hl.is_defined(ht_snv_age_accum.HL) & (ht_snv_age_accum.HL < 0.95) & (ht_snv_age_accum.HL >= 0.01) & (ht_snv_age_accum.AD[1] > ht_snv_age_accum.pois_cutoff), hl.agg.sum(1 / ht_snv_age_accum.HL)),
                                                            snv_ADpois95_1minHL = hl.agg.filter(hl.is_defined(ht_snv_age_accum.HL) & (ht_snv_age_accum.HL < 0.95) & (ht_snv_age_accum.HL >= 0.01) & (ht_snv_age_accum.AD[1] > ht_snv_age_accum.pois_cutoff), hl.agg.sum(1 - ht_snv_age_accum.HL))).select_globals()
        ht_snv_age_accum_count = ht_snv_age_accum_count.annotate(snv_meanHL = hl.if_else(hl.is_nan(ht_snv_age_accum_count.snv_meanHL), 0, ht_snv_age_accum_count.snv_meanHL))
        ht_snv_age_accum_count = ht_snv_age_accum_count.checkpoint(get_final_munged_snvcount_age_accum_path(version, HL), overwrite=True)
        ht_snv_age_accum_count.export(get_final_munged_snvcount_age_accum_path(version, HL, 'tsv'))
    
    return ht_snv_age_accum_count



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