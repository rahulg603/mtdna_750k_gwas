#!/usr/bin/env python3
import hail as hl
import pandas as pd
import subprocess
import os, re

from paths import *


def get_final_annotated_variant_path():
    return f'{BUCKET}/final_v7_merged_callset/all_v7_callset_20241030/vcf/filt_annotated/annotated_combined_processed_flat.tsv.bgz'


def get_final_munged_case_only_hl_path(num_to_keep, type):
    base_path = f'{BUCKET}/final_v7_merged_callset/all_v7_callset_20241030/munged/'
    if type == 'ht':
        return base_path + f'case_only_calls_low_hl_fullqc_N_{str(num_to_keep)}.ht'
    elif type == 'tsv':
        return base_path + f'case_only_calls_low_hl_fullqc_long_format_N_{str(num_to_keep)}.tsv'


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


def get_case_only_mtdna_callset(num_to_keep=310, overwrite=False):
    """ Here we implement the munging pipeline in Python and Pandas.
    """
    if hl.hadoop_exists(f"{get_final_munged_case_only_hl_path(num_to_keep, 'ht')}/_SUCCESS") and not overwrite:
        ht = hl.read_table(get_final_munged_case_only_hl_path(num_to_keep, 'ht'))
    
    else:
        df = pd.read_csv(get_final_annotated_variant_path(), sep='\t', 
                         compression='gzip',
                         dtype={'locus': 'str', 'alleles': 'str', 's': 'str',
                                'rsid': 'str', 'variant':'str', 'batch':'str',
                                'common_low_heteroplasmy': 'boolean',
                                'hap_defining_variant': 'boolean',
                                'region': 'str',
                                'variant_context': 'str',
                                'dp_mean': 'float64',
                                'mq_mean': 'float64', 'tlod_mean': 'float64',
                                'AF_hom': 'float64', 'AF_het': 'float64',
                                'AC_hom': 'int64', 'AC_het': 'int64',
                                'max_hl': 'float64', 'dp_mean': 'float64',
                                'major_haplogroup':'str', 'hap':'str', 
                                'wgs_median_coverage':'float64', 'mt_mean_coverage':'float64',
                                'mito_cn':'float64', 'age':'float64', 'pop':'str',
                                'DP':'float64', 'AD_ref':'float64','AD_alt':'float64',
                                'MQ':'float64','TLOD':'float64', 'GT':'str',
                                'OriginalSelfRefAlleles':'str', 'SwappedFieldIDs':'str', 'FT':'str', 'FT_LIFT':'str',
                                'artifact_prone':'boolean', 'lifted':'boolean', 'fail_gt':'boolean', 'missing_call':'boolean'})
        df_for_enumeration = df[(~df['HL'].isna()) & (df['HL'] < 0.95) & (df['common_low_heteroplasmy'])]
        counted_variants = df_for_enumeration.groupby('variant')['HL'].count().sort_values(ascending=False)
        variants_to_extract = list(counted_variants[counted_variants > num_to_keep].index)

        variants_to_analyze = df[df['variant'].isin(variants_to_extract) & (~df['HL'].isna()) & (df['HL'] < 0.95)]
        variants_to_analyze.to_csv(get_final_munged_case_only_hl_path(num_to_keep, 'tsv'), sep='\t', index=False)
        
        pivoted_vars = variants_to_analyze.pivot(index='s', columns='variant',values='HL').reset_index()
        alls = pd.DataFrame({'s': list(set(df['s']))})
        data_HL = pd.merge(pivoted_vars, alls, on='s', how='right')
        backbone = data_HL[['s']].copy()
        for vari in variants_to_extract:
            this_name = re.sub('[:,]','_', vari)
            backbone.loc[:,this_name] = data_HL[vari]

        backbone.to_csv(f'{TEMP_PATH}pivoted_flat_file_lowHLqc.tsv', sep='\t', index=False)

        ht = hl.import_table(f'{TEMP_PATH}pivoted_flat_file_lowHLqc.tsv', impute=True, missing='')
        ht = ht.annotate(s = hl.str(ht.s)).key_by('s')
        ht = ht.repartition(50).checkpoint(get_final_munged_case_only_hl_path(num_to_keep, 'ht'))
    
    return ht


def get_snv_count_phenotype():
    ht = hl.import_table(get_final_annotated_variant_path(),
                         types={'locus': hl.tstr, 'alleles': hl.tstr, 's': hl.tstr,
                                'rsid': hl.tstr, 'variant':hl.tstr, 'batch':hl.tstr,
                                'common_low_heteroplasmy': hl.tbool,
                                'hap_defining_variant': hl.tbool,
                                'region': hl.tstr,
                                'variant_context': hl.tstr,
                                'dp_mean': hl.tfloat64,
                                'mq_mean': hl.tfloat64, 'tlod_mean': hl.tfloat64,
                                'AF_hom': hl.tfloat64, 'AF_het': hl.tfloat64,
                                'AC_hom': hl.tint64, 'AC_het': hl.tint64,
                                'max_hl': hl.tfloat64, 'dp_mean': hl.tfloat64,
                                'major_haplogroup':hl.tstr, 'hap':hl.tstr, 
                                'wgs_median_coverage':hl.tfloat64, 'mt_mean_coverage':hl.tfloat64,
                                'mito_cn':hl.tfloat64, 'age':hl.tfloat64, 'pop':hl.tstr,
                                'DP':hl.tfloat64, 'AD_ref':hl.tfloat64,'AD_alt':hl.tfloat64,
                                'MQ':hl.tfloat64,'TLOD':hl.tfloat64, 'GT':hl.tstr, 'HL':hl.tfloat64,
                                'OriginalSelfRefAlleles':hl.tstr, 'SwappedFieldIDs':hl.tstr, 'FT':hl.tstr, 'FT_LIFT':hl.tstr,
                                'artifact_prone':hl.tbool, 'lifted':hl.tbool, 'fail_gt':hl.tbool, 'missing_call':hl.tbool}, min_partitions=50)
    all_s_table = ht.select('s').key_by('s').distinct()
    ht_heteroplasmies = ht.filter(hl.is_defined(ht.HL) & (ht.HL < 0.95)).repartition(40)
    ht_heteroplasmies = ht_heteroplasmies.annotate(alleles = ht_heteroplasmies.alleles.split(','))
    ht_heteroplasmies = ht_heteroplasmies.filter(~hl.is_indel(ht_heteroplasmies.alleles[0], ht_heteroplasmies.alleles[1]))
    ht_snv_count = ht_heteroplasmies.group_by(ht_heteroplasmies.s).aggregate(N = hl.agg.count())
    s_not_found = all_s_table.anti_join(ht_snv_count).annotate(N = hl.int64(hl.literal(0)))
    ht_snv_count = hl.Table.union(ht_snv_count, s_not_found)
    return ht_snv_count.rename({'N':'snv_count_qcpass'})
