#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GxG detection heart_rate_tmp
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Author: Saul Pierotti
Mail: saul@ebi.ac.uk
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

process read_pheno_covar {
    label "r_tidyverse_datatable"

    input:
        path pheno
        path covar

    output:
        path "pheno_covar.csv.gz"

    script:
        """
        #!/usr/bin/env Rscript

        library("data.table")
        library("tidyverse")

        pheno <- fread("${pheno}")
        covar <- fread("${covar}")

        pheno <- pheno[, .(individual = IID, heart_rate_avg_21C, heart_rate_avg_28C, heart_rate_avg_35C)]
        pheno <- melt(
            pheno,
            id.vars = "individual",
            measure.vars = c("heart_rate_avg_21C", "heart_rate_avg_28C", "heart_rate_avg_35C"),
            value.name = "heart_rate",
            variable.name = "temperature"
        )
        pheno[, temperature := str_remove(temperature, "heart_rate_avg_")]
        # the kronoecker product creates id of the kind individual:temperature
        pheno[, full_id := sprintf("%s:%s", individual, temperature)]

        covar <- covar[, .(individual = IID, phenotyping_plate_id, cross_id)]

        df <- merge(pheno, covar, by = "individual")
        fwrite(df, "pheno_covar.csv.gz")
        """
}

process get_qtl_tests_and_models {
    label "r_tidyverse_datatable"

    input:
        path qtls
        path formulas
        path testing_scheme

    output:
        path "qtl_pairs.csv.gz", emit: pairs
        path "qtl_tests.csv.gz", emit: tests
        path "qtl_models.csv.gz", emit: models

    script:
        """
        #!/usr/bin/env Rscript

        library("data.table")

        qtls <- fread("${qtls}")
        formulas <- readRDS("${formulas}")
        testing_scheme <- readRDS("${testing_scheme}") |> as.data.table()

        qtls <- qtls[
            qtl_type == "temp" & well_behaved == "yes",
            .(locus_id1 = locus_id, chr1 = chr, lead_snp_id1 = lead_snp_id)
        ]
        # all possible n choose 2 combinations of 2 loci
        all_comb <- combn(1:nrow(qtls), 2)
        qtls <- cbind(
            qtls[all_comb[1,]],
            qtls[all_comb[2,], .(locus_id2 = locus_id1, lead_snp_id2 = lead_snp_id1, chr2 = chr1)]
        )

        qtl_models <- merge(
            qtls[, k := ""],
            data.table(model = names(formulas), k = ""),
            allow.cartesian = TRUE,
            by = "k"
        )
        

        qtl_tests <- merge(
            qtls[, k := ""],
            testing_scheme[, k := ""],
            allow.cartesian = TRUE,
            by = "k"
        )[, k := NULL][]

        qtls[, k := NULL]
        qtl_models[, k := NULL]
        qtl_tests[, k := NULL]

        fwrite(qtls, "qtl_pairs.csv.gz")
        fwrite(qtl_tests, "qtl_tests.csv.gz")
        fwrite(qtl_models, "qtl_models.csv.gz")
        """
}

process get_formulas_and_testing_scheme {
    label "r_tidyverse_datatable"

    output:
        path "formulas.rds", emit: formulas
        path "testing_scheme.rds", emit: testing_scheme

    script:
        """
        #!/usr/bin/env Rscript

        formulas <- list(
            gxgxe_gxg_gxe = formula(heart_rate ~ 1 + cross_id*temperature + phenotyping_plate_id + temperature + snp1 + snp1:temperature + snp2 + snp2:temperature + snp1:snp2 + snp1:snp2:temperature),
            gxg_gxe = formula(heart_rate ~ 1 + cross_id*temperature + phenotyping_plate_id + temperature + snp1 + snp1:temperature + snp2 + snp2:temperature + snp1:snp2),
            gxe = formula(heart_rate ~ 1 + cross_id*temperature + phenotyping_plate_id + temperature + snp1 + snp2 + snp1:temperature + snp2:temperature)
        )

        testing_scheme <- testing_scheme <- list(
            model         = c("gxg_gxe", "gxgxe_gxg_gxe"),
            reduced_model = c("gxe", "gxg_gxe")
        )

        saveRDS(formulas, "formulas.rds")
        saveRDS(testing_scheme, "testing_scheme.rds")
        """
}

process make_freq {
    label "plink2"

    input:
        path vcf

    output:
        path "freq.afreq.zst"

    script:
        """
        plink2 \\
            --threads ${task.cpus} \\
            --memory ${task.memory.getMega()} \\
            --vcf ${vcf} \\
            --min-alleles 2 \\
            --max-alleles 2 \\
            --set-missing-var-ids @_#_\\\$r_\\\$a \\
            --new-id-max-allele-len 10 missing \\
            --out freq \\
            --chr-set ${params.n_chr} \\
            --freq zs
        """
}

process make_pgen {
    label "plink2"

    input:
        path vcf

    output:
        tuple(
            path("pgen.pgen"),
            path("pgen.psam"),
            path("pgen.pvar.zst")
        )

    script:
        """
        plink2 \\
            --threads ${task.cpus} \\
            --memory ${task.memory.getMega()} \\
            --set-all-var-ids @_#_\\\$r_\\\$a \\
            --min-alleles 2 \\
            --max-alleles 2 \\
            --out pgen \\
            --chr-set ${params.n_chr} \\
            --vcf ${vcf} dosage=DS \\
            --make-pgen vzs fill-missing-from-dosage erase-dosage
        """
}

process make_grm {
    label "plink2"
    tag "${meta.id}"

    input:
        tuple(
            val(meta),
            path(pgen),
            path(psam),
            path(pvar),
            path(freq)
        )

    output:
        tuple(
            val(meta),
            path("grm_loco_${meta.id}.rel.id"),
            path("grm_loco_${meta.id}.rel.bin")
        )

    script:
        """
        plink2 \
            --pgen $pgen \\
            --psam $psam \\
            --pvar $pvar \\
            --threads ${task.cpus} \\
            --memory ${task.memory.getMega()} \\
            --not-chr "${meta.chr1} ${meta.chr2}" \
            --read-freq $freq \
            --out "grm_loco_${meta.id}" \
            --maf 0.01 \
            --make-rel bin
        """
}

process get_qtl_matrices {
    label "r_tidyverse_datatable"
    tag "${meta.id}"

    input:
        tuple(
            val(meta),
            path(grm_id),
            path(grm_bin),
            path(pheno_covar),
            path(pgen),
            path(psam),
            path(pvar),
            path(formulas)
        )

    output:
        tuple(
            val(meta),
            path("${meta.id}.qtl_matrices.rds")
        )

    script:
        """
        #!/usr/bin/env Rscript

        library("data.table")

        get_snp <- function(var_id){
            pvar <- pgenlibr::NewPvar("${pvar}")
            pgen <- pgenlibr::NewPgen("${pgen}", pvar = pvar)
            ret <- fread("${psam}")[, .(individual = `#IID`)]
            ret[["snp"]] <- pgenlibr::Buf(pgen)
            var_num <- pgenlibr::GetVariantsById(pvar, var_id)
            pgenlibr::Read(pgen, ret[["snp"]], var_num, allele_num = 2L)
            ret[["dominance"]] = (ret[["snp"]] == 1)
            return(ret)
        }

        get_design_matrix <- function(the_formula, model_frame){
            X <- model.matrix(the_formula, data = model_frame)
            return(X)
        }

        read_K <- function(grm_id, grm_bin){
            samples_K <- read.table(grm_id, header = FALSE, check.names = FALSE)[,1]
            K <- matrix(
                readBin(grm_bin, what = "numeric", n = length(samples_K) ** 2),
                ncol = length(samples_K)
            )
            colnames(K) <- samples_K
            rownames(K) <- samples_K
            stopifnot(sum(is.na(K)) == 0)
            return(K)
        }

        get_K_ind <- function(df, model_frame){
            # expand the factors in a dummy encoded matrix with n_samples rows and n_groups columns
            # needs to be df and not model_frame because the latter does not contain indivudual
            Z <- model.matrix(~ 0 + individual, data = df)
            # get a block diagonal square matrix with n_samples rows and columns and
            # 1 for samples in the same group, 0 for samples in different groups
            K <- Z %*% t(Z)
            
            match_v <- match(rownames(model_frame), rownames(K))
            K <- K[match_v, match_v]
            
            return(K)
        }

        get_K_list <- function(grm_id, grm_bin, model_frame, df){
            K_grm <- read_K(grm_id, grm_bin)
            possible_temps <- model_frame[["temperature"]] |> unique()
            n_temps <- length(possible_temps)
            # reletedness among environments is diagonal
            Kt_indep <- diag(n_temps)
            rownames(Kt_indep) <- possible_temps
            colnames(Kt_indep) <- possible_temps
            # relatedness among environment is full
            Kt_full <- matrix(1, n_temps, n_temps)
            rownames(Kt_full) <- possible_temps
            colnames(Kt_full) <- possible_temps
            # kronecker products of the relationships among environments and genetic relatedness
            # order does not matter since I match samples afterwards (matters for the row and colnames)
            K_grm_indep <- kronecker(K_grm, Kt_indep, make.dimnames = TRUE) # this is used to model s2_gxe
            K_grm_full <- kronecker(K_grm, Kt_full, make.dimnames = TRUE) # this is used to model s2_g
            # s2_g and s2_gxe together create the "compound symmetry" model for GxE where an overall variance s2_g is estimated and a within-environment variance s2_gxe
            # more complex models are possible but probably overkill
            match_v <- match(rownames(model_frame), rownames(K_grm_indep))
            # one match is sufficient, Kt_indep and Kt_full have the same names
            ret <- list(
                K_grm_indep = K_grm_indep[match_v, match_v],
                K_grm_full = K_grm_full[match_v, match_v],
                K_ind = get_K_ind(df, model_frame) # correlation among sample duplicates
            )
            return(ret)
        }

        pheno_covar <- fread("${pheno_covar}")
        formulas <- readRDS("${formulas}")
        the_formula <- formulas[["${meta.model}"]]
        snp1 <- get_snp("${meta.lead_snp_id1}")
        snp2 <- get_snp("${meta.lead_snp_id2}")
        df <- merge(pheno_covar, snp1[, .(snp1 = snp, individual)], by = "individual")
        df <- merge(df, snp2[, .(snp2 = snp, individual)], by = "individual")
        df <- as.data.frame(df)
        rownames(df) <- df[["full_id"]]
        model_frame <- model.frame(the_formula, data = df)
        K_list <- get_K_list("${grm_id}", "${grm_bin}", model_frame, df)
        y <- model.response(model_frame)
        X <- get_design_matrix(the_formula, model_frame)
        res <- list(
            y = y,
            X = X,
            K_list = K_list
        )

        saveRDS(res, "${meta.id}.qtl_matrices.rds")
        """
}

process fit_mixed_model {
    label "r_gaston"
    tag "${meta.id}"

    input:
        tuple(
            val(meta),
            path(qtl_matrices)
        )

    output:
        path "${meta.id}.mm_fit.rds"

    script:
        """
        #!/usr/bin/env Rscript

        qtl_matrices <- readRDS("${qtl_matrices}")
        X <- qtl_matrices["X"]]
        y <- qtl_matrices["y"]]
        K_list <- qtl_matrices[["K_list"]]
        fit <- gaston::lmm.aireml(
            Y = y,
            X = X,
            K = K_list,
            verbose = TRUE
        )

        saveRDS(fit, "${meta.id}.mm_fit.rds")
        """
}


workflow {
    read_pheno_covar ( params.pheno, params.covar )
    get_formulas_and_testing_scheme ()
    get_formulas_and_testing_scheme.out.formulas.set { formulas }
    get_formulas_and_testing_scheme.out.testing_scheme.set { testing_scheme }
    get_qtl_tests_and_models ( params.qtls, formulas, testing_scheme )
    get_qtl_tests_and_models.out.pairs
        .splitCsv ( header: true )
        .map {
            it.id = it.locus_id1 + "_" + it.locus_id2
            return ( it )
        }
        .set { qtl_pairs }
    get_qtl_tests_and_models.out.models
        .splitCsv ( header: true )
        .map {
            it.id = it.locus_id1 + "_" + it.locus_id2 + "_" + it.model
            return ( it )
        }
        .set { qtl_models }
    make_freq ( params.freq )
    make_pgen ( params.vcf )
    qtl_pairs.combine ( make_pgen.out ).combine ( make_freq.out ).set { make_grm_in_ch }
    make_grm ( make_grm_in_ch )
    make_grm.out
        .map { meta, id, bin -> [[meta.locus_id1, meta.locus_id2], meta, id, bin] }
        .combine ( qtl_models.map { meta -> [[meta.locus_id1, meta.locus_id2], meta] }, by: 0 )
        .map { match_tuple, meta1, id, bin, meta2 -> [meta2, id, bin] }
        .combine ( read_pheno_covar.out )
        .combine ( make_pgen.out )
        .combine ( formulas )
        .set { get_qtl_matrices_in_ch }
    get_qtl_matrices ( get_qtl_matrices_in_ch )
    get_qtl_matrices.out
        // fit variance components only on the background model
        .filter { meta, qtl_mat -> meta.model == "gxe" }
        .set { fit_mixed_model_in_ch }
    fit_mixed_model ( fit_mixed_model_in_ch )
}
