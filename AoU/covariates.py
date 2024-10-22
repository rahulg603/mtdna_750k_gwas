#!/usr/bin/env python3
import hail as hl
import pandas as pd
import os
import re

from AoU.paths import *
from utils.SaigeImporters import *


def generate_indicator(ht, col, baseline_item):
    assert(ht[col].dtype == hl.tstr)

    all_values = ht.aggregate(hl.agg.collect_as_set(ht[col]))
    if baseline_item not in all_values:
        raise ValueError('Item to be used as baseline is not present in data.')
    return ht.annotate(**{f'{col}_{x}': ht[col] == x for x in all_values if x != baseline_item})


def load_ancestry_data(use_drc_data=True):
    # if use_drc_data is True, will return the first 15 covariates as produced by DRC.
    # if use_drc_data is False, will return the first 20 covariates as produced for AllxAllofUs
    if use_drc_data:
        npcs = 15
        ancestry_pred = hl.import_table(get_ancestry_flat_file(),
                                        key="research_id", 
                                        impute=True, 
                                        types={"research_id":"tstr","pca_features":hl.tarray(hl.tfloat)},
                                        min_partitions=50)
        ancestry_pred = ancestry_pred.select(pop = ancestry_pred.ancestry_pred,
                                             **{f'PC{str(idx+1)}': ancestry_pred.pca_features[idx] for idx in range(0,npcs+1)})
        ancestry_pred = ancestry_pred.rename({'research_id': 's'}).key_by('s')
        return npcs, ancestry_pred
    else:
        # TODO AxAoU data
        raise NotImplementedError('AxAoU ancestry data not implemented yet.')


def load_aou_genotyping_qc_data():
    ht_qc_flat = hl.import_table(get_genomic_metric_flat_file(),
                                 impute=True, 
                                 types={"research_id":"tstr",
                                        "verify_bam_id2_contamination":"tfloat"},
                                 min_partitions=50,
                                 missing=['NA',''])
    ht_qc_flat = ht_qc_flat.annotate(isFemale = hl.if_else(ht_qc_flat.sex_at_birth == 'F', 1, 
                                                hl.if_else(ht_qc_flat.sex_at_birth == 'M', 0, hl.missing(hl.tint32))))
    ht_qc_flat = ht_qc_flat.rename({'research_id': 's'}).key_by('s')
    return ht_qc_flat


def load_axaou_sample_qc():
    ht = hl.import_table(get_all_by_all_qc_samples(),
                         impute=False,
                         types={"person_id":"tstr"})
    ht = ht.rename({'person_id': 's'}).key_by('s')
    return ht


def load_aou_flagged_samples():
    ht = hl.import_table(get_aou_flagged_samples(),
                         impute=True,
                         types={"s":"tstr",
                                "probabilities":hl.tarray(hl.tfloat),
                                "pca_features":hl.tarray(hl.tfloat),
                                "qc_metrics_filters":hl.tarray(hl.tstr)},
                        key='s')
    return ht


def get_all_demographics(overwrite=False, use_drc_ancestry_data=True):
    covar_path = get_demographics_path(COVAR_PATH, use_drc_ancestry_data)
    if not overwrite and hl.hadoop_exists(os.path.join(covar_path, '_SUCCESS')):
        sample_covariates = hl.read_table(covar_path)
    else:
        npcs, ancestry_pred = load_ancestry_data(use_drc_ancestry_data)
        dataset = os.getenv('WORKSPACE_CDR')
        person_sql = f"""
        SELECT  person.person_id,
                person.birth_datetime,
                p_location_concept.concept_name as loc,
                p_site_concept.concept_name as site
            FROM
                `{dataset}.person` person 
            LEFT JOIN
                `{dataset}.concept` p_location_concept 
                    on person.location_id = p_location_concept.CONCEPT_ID 
            LEFT JOIN
                `{dataset}.concept` p_site_concept 
                    on person.care_site_id = p_site_concept.CONCEPT_ID
            WHERE
                person.PERSON_ID IN (
                    select
                        person_id  
                    from
                        `{dataset}.cb_search_person` cb_search_person  
                    where
                        cb_search_person.person_id in (
                            select
                                person_id 
                            from
                                `{dataset}.cb_search_person` p 
                            where
                                has_whole_genome_variant = 1 
                        ) 
                    )"""

        wgs_demog = pd.read_gbq(person_sql, dialect="standard")
        wgs_demog['birth_year'] = pd.DatetimeIndex(wgs_demog.birth_datetime).year
        age_ht = hl.Table.from_pandas(wgs_demog[['person_id', 'birth_year']])
        age_ht = age_ht.annotate(s = hl.str(age_ht.person_id)).key_by('s')

        ht_qc_flat = load_aou_genotyping_qc_data()
        axaou_sample_qc = load_axaou_sample_qc()
        flagged_samples = load_aou_flagged_samples()

        # munge covariates
        sample_covariates = ht_qc_flat.select(isFemale = ht_qc_flat.isFemale, 
                                              nucdna_mean_coverage = ht_qc_flat.mean_coverage,
                                              site_id = ht_qc_flat.site_id,
                                              sample_collection_date = ht_qc_flat.biosample_collection_date,
                                              sample_source = ht_qc_flat.sample_source,
                                              ploidy = ht_qc_flat.dragen_sex_ploidy,
                                              contamination = hl.struct(dragen = ht_qc_flat.dragen_contamination,
                                                                     verify_bam_id2 = ht_qc_flat.verify_bam_id2_contamination))
        sample_covariates = sample_covariates.annotate(sample_collection_split = sample_covariates.sample_collection_date.split('-'))
        sample_covariates = sample_covariates.annotate(sample_collection_year = hl.if_else(hl.len(sample_covariates.sample_collection_split) == 3,
                                                                                           hl.int32(sample_covariates.sample_collection_split[0]),
                                                                                           hl.missing('int32'))).drop('sample_collection_split')
        
        # add birth year and compute age
        sample_covariates = sample_covariates.annotate(birth_year = age_ht[sample_covariates.s].birth_year)
        sample_covariates = sample_covariates.annotate(age = sample_covariates.sample_collection_year - sample_covariates.birth_year)

        # add GWAS covariates
        sample_covariates = sample_covariates.annotate(age_isFemale = sample_covariates.age * sample_covariates.isFemale,
                                                       age2 = sample_covariates.age**2)
        sample_covariates = sample_covariates.annotate(age2_isFemale = sample_covariates.age2 * sample_covariates.isFemale)
        sample_covariates = sample_covariates.annotate(ancestry = hl.struct(pop = ancestry_pred[sample_covariates.s].pop,
                                                                            **{f'PC{str(idx+1)}': ancestry_pred[sample_covariates.s][f'PC{str(idx+1)}'] for idx in range(0,npcs+1)}))

        # add flag for passing axaou
        sample_covariates = sample_covariates.annotate(pass_axaou_qc = hl.is_defined(axaou_sample_qc[sample_covariates.s]))

        # add flag from AoU DRC
        sample_covariates = sample_covariates.annotate(pass_aou_drc_qc = ~hl.is_defined(flagged_samples[sample_covariates.s]),
                                                       failed_flags_drc_qc = flagged_samples[sample_covariates.s].qc_metrics_filters)
        sample_covariates = sample_covariates.annotate(pass_qc = sample_covariates.pass_axaou_qc & sample_covariates.pass_aou_drc_qc)

        # output
        sample_covariates = sample_covariates.repartition(100).checkpoint(covar_path, overwrite=True)
    
    return sample_covariates


def get_gwas_covariates(overwrite=False, use_drc_ancestry_data=True):
    baseline_covar_path = get_base_covariates_path(COVAR_PATH, drc=use_drc_ancestry_data)
    if overwrite or not hl.hadoop_exists(os.path.join(baseline_covar_path, '_SUCCESS')):
        # now producing baseline covariates, which includes making indicator variables for sequencing center and keeping PCs 1-16, age, age2, sex, age2_sex
        sample_covariates = get_all_demographics(overwrite=overwrite, use_drc_ancestry_data=use_drc_ancestry_data)
        baseline_covar = generate_indicator(sample_covariates, col='site_id', baseline_item='bi').key_by('s')
        dict_update = {}
        dict_update.update({x: baseline_covar['ancestry'][x] for x in baseline_covar.ancestry.keys() if re.search('^PC[0-9]{1,2}$', x)})
        dict_update.update({x: baseline_covar[x] for x in baseline_covar.row.keys() if re.search('^site_id_.+$', x)})
        dict_update.update({'sex': baseline_covar['isFemale'],
                                 'age': baseline_covar['age'],
                                 'age_sex': baseline_covar['age_isFemale'],
                                 'age2': baseline_covar['age2'],
                                 'age2_sex': baseline_covar['age2_isFemale'],
                                 'pop': baseline_covar['ancestry']['pop']})
        baseline_covar = baseline_covar.select(**dict_update)
        baseline_covar = baseline_covar.repartition(50).checkpoint(baseline_covar_path, overwrite=True)
    else:
        baseline_covar = hl.read_table(baseline_covar_path)

    return baseline_covar