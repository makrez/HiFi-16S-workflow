nextflow.enable.dsl = 2

include {
    taxonomy_nb_assign
    taxonomy_nb_best
    taxonomy_nb_merge
} from '../modules/taxonomy_nb'

include {
    taxonomy_add_md5 as taxonomy_nb_add_md5
    taxonomy_summary as taxonomy_nb_summary
} from '../modules/taxonomy_common'


workflow TAXONOMY_NB {

    take:
    nb_inputs_ch
    asv_table_tsv
    db_priority

    main:

    taxonomy_nb_assign(nb_inputs_ch)

    taxonomy_nb_best(
        taxonomy_nb_assign.out.nb_tax
            .map { db_name, taxonomy -> taxonomy }
            .collect(),
        db_priority
    )

    taxonomy_nb_merge(
        taxonomy_nb_best.out.best_nb_tax,
        asv_table_tsv
    )

    taxonomy_nb_add_md5(
        taxonomy_nb_merge.out.merged_no_id,
        'best_nb'
    )

    taxonomy_nb_summary(
        taxonomy_nb_add_md5.out.tax_with_md5,
        'best_nb'
    )

    emit:
    nb_tax                  = taxonomy_nb_assign.out.nb_tax
    best_nb_tax             = taxonomy_nb_best.out.best_nb_tax
    best_nb_tax_with_db     = taxonomy_nb_best.out.best_nb_tax_with_db
    final_nb_table          = taxonomy_nb_add_md5.out.tax_with_md5
    assignment_summary      = taxonomy_nb_summary.out.assignment_summary
    deepest_rank_summary    = taxonomy_nb_summary.out.deepest_rank_summary
}
