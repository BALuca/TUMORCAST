#!/usr/bin/env Rscript
# Extend the dataset59 inventory: actual local drug-cell support, not inferred
# Cartesian products of drug and cell catalogs. Earlier outputs are preserved.
suppressPackageStartupMessages({library(data.table); library(qs)})
setDTthreads(2L)
root <- '/work/results/program_pilot/59_mechanistic_evidence_datasets'
out <- file.path(root, 'data')
drug_key <- function(x) tolower(gsub('[^a-zA-Z0-9]', '', x))
cell_key <- function(x) sub('^NCI', '', toupper(gsub('[^a-zA-Z0-9]', '', sub('_.*', '', x))))
aliases <- list(tamoxifen = c('tamoxifen','tamoxifen citrate'), mefloquine = c('mefloquine','mefloquine hydrochloride'),
  withaferin = c('withaferin a'), trametinib = c('trametinib','trametinib DMSO TF solvate','GSK1120212'),
  lestaurtinib = c('lestaurtinib','CEP701'), midostaurin = c('midostaurin','PKC412'),
  sotorasib = c('sotorasib','AMG510'), selumetinib = c('selumetinib','AZD6244'),
  pictilisib = c('pictilisib','GDC0941'), jnkin8 = c('JNKIN8','JNKinVIII'),
  lapatinib = c('lapatinib','lapatinib ditosylate','lapatinib ditosylate monohydrate','GW572016'))
alias_table <- rbindlist(lapply(names(aliases), function(nm) data.table(drug = nm, alias = aliases[[nm]])))
alias_table[, key := drug_key(alias)]
stopifnot(!anyDuplicated(alias_table$key))
canonical <- function(x) alias_table$drug[match(drug_key(x), alias_table$key)]
targets <- rbindlist(list(
  CJ(study = 'Diaz2020', cell = c('MCF7','LNCAP'), drug = c('tamoxifen','mefloquine','withaferin')),
  CJ(study = 'Macaya2023_molecular', cell = 'H1792', drug = c('trametinib','lestaurtinib','midostaurin')),
  data.table(study = 'Macaya2023_genetic_control', cell = c('H2009','H358','H358'), drug = c('trametinib','trametinib','sotorasib')),
  CJ(study = 'vanBentum2025', cell = c('HCT116','DLD1','CACO2'), drug = c('selumetinib','pictilisib','jnkin8')),
  data.table(study = 'ERBB_rescue', cell = 'SKBR3', drug = 'lapatinib')))
targets[, cell := cell_key(cell)]

# Prepared Tahoe profiles retain actual drug-cell availability. The catalog
# separately establishes whether a parent drug exists in this pinned PDEX set.
tahoe_catalog_path <- '/work/cache/tahoe_prism_full_chain/tahoe_pdex_selected/perturbations.tsv'
tahoe_catalog <- fread(tahoe_catalog_path)
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
panel <- qread(panel_path)
tahoe <- copy(panel$metadata)
tahoe[, `:=`(source_drug = drug, drug = canonical(drug), cell = cell_key(recipient))]
tahoe <- tahoe[!is.na(drug)]

# PRISM requires a finite response in the actual cell, not just a treatment-info
# entry. Keep the discrete measured doses; do not interpolate for this audit.
prism_dir <- '/work/cache/tahoe_prism_full_chain/prism'
prism_treatment <- fread(file.path(prism_dir, 'secondary-screen-replicate-collapsed-treatment-info.csv'))
prism_treatment[, drug := canonical(name)]
prism_treatment <- prism_treatment[!is.na(drug)]
prism_cells <- fread(file.path(prism_dir, 'secondary-screen-cell-line-info.csv'))
prism_cells[, cell := cell_key(ccle_name)]
prism_path <- file.path(prism_dir, 'secondary-screen-replicate-collapsed-logfold-change.csv')
first <- names(fread(prism_path, nrows = 0))[1]
prism_wide <- fread(prism_path, select = c(first, prism_treatment$column_name))
setnames(prism_wide, first, 'depmap_id')
prism_wide <- prism_wide[depmap_id %chin% prism_cells[cell %chin% targets$cell, depmap_id]]
prism <- melt(prism_wide, id.vars = 'depmap_id', variable.name = 'column_name', value.name = 'response')
prism <- merge(prism, prism_treatment[, .(column_name, drug, source_drug = name, dose)], by = 'column_name')
prism <- merge(prism, unique(prism_cells[, .(depmap_id, cell)]), by = 'depmap_id')[is.finite(response)]

lincs_path <- '/work/v20_pilot_design_2026-10-09/source/GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz'
lincs <- fread(lincs_path)
lincs[, `:=`(drug = canonical(pert_iname), cell = cell_key(cell_id))]
lincs <- lincs[!is.na(drug) & pert_type == 'trt_cp']

# Original and validation Jaaks matrices contain observed combinations. Preserve
# both actual partners, rather than treating each listed drug as all pairings.
jaaks_paths <- file.path('/work/jaaks_2022/raw', c('Original screen_All tissues_fitted data.csv',
  'Validation screen_All tissues_fitted data.csv'))
jaaks <- rbindlist(lapply(jaaks_paths, function(f) {
  x <- fread(f, select = c('CELL_LINE_NAME','ANCHOR_NAME','LIBRARY_NAME','ANCHOR_CONC','LIBRARY_CONC','Synergy'))
  x[, source := basename(f)]
  x
}))
jaaks[, `:=`(cell = cell_key(CELL_LINE_NAME), drug_A = canonical(ANCHOR_NAME), drug_B = canonical(LIBRARY_NAME))]
jaaks <- jaaks[!is.na(drug_A) | !is.na(drug_B)]

almanac_names_path <- '/work/almanac_conditionals/raw/ComboCompoundNames_small.txt'
almanac_names <- fread(almanac_names_path, col.names = c('nsc','source_drug'))
almanac_names[, drug := canonical(source_drug)]
almanac_path <- '/work/almanac_conditionals/raw/ComboDrugGrowth_Nov2017.csv'
almanac <- fread(almanac_path, select = c('CELLNAME','NSC1','NSC2','CONC1','CONC2','VALID','PERCENTGROWTH'))
almanac[, `:=`(cell = cell_key(CELLNAME), drug_A = almanac_names$drug[match(NSC1, almanac_names$nsc)],
  drug_B = almanac_names$drug[match(NSC2, almanac_names$nsc)])]
almanac <- almanac[(!is.na(drug_A) | !is.na(drug_B)) & is.finite(PERCENTGROWTH)]

rows <- list()
for (i in seq_len(nrow(targets))) {
  z <- targets[i]
  th <- tahoe[cell == z$cell & drug == z$drug]
  pr <- prism[cell == z$cell & drug == z$drug]
  li <- lincs[cell == z$cell & drug == z$drug]
  ja <- jaaks[cell == z$cell & (drug_A == z$drug | drug_B == z$drug)]
  al <- almanac[cell == z$cell & (drug_A == z$drug | drug_B == z$drug)]
  rows[[i]] <- cbind(z, data.table(
    Tahoe_drug_in_PDEX_catalog = z$drug %in% canonical(tahoe_catalog$drug),
    Tahoe_prepared_profiles = nrow(th), Tahoe_doses_uM = paste(sort(unique(th$dose)), collapse = ';'),
    PRISM_finite_measurements = nrow(pr), PRISM_doses_uM = paste(sort(unique(pr$dose)), collapse = ';'),
    LINCS_GSE70138_samples = nrow(li), LINCS_doses = paste(sort(unique(li$pert_dose)), collapse = ';'),
    LINCS_times_h = paste(sort(unique(li$pert_time)), collapse = ';'),
    Jaaks_combination_fit_rows_involving_drug = nrow(ja),
    ALMANAC_rows_involving_drug = nrow(al)))
}
coverage <- rbindlist(rows)
fwrite(coverage, file.path(out, 'resource_drug_cell_overlap.tsv'), sep = '\t')
fwrite(alias_table, file.path(out, 'resource_overlap_aliases.tsv'), sep = '\t')
fwrite(unique(lincs[cell %chin% targets$cell, .(cell, drug, source_drug = pert_iname, dose = pert_dose,
  dose_unit = pert_dose_unit, time_h = pert_time, inst_id)]), file.path(out, 'resource_LINCS_matching_samples.tsv'), sep = '\t')
fwrite(prism, file.path(out, 'resource_PRISM_matching_measurements.tsv'), sep = '\t')
fwrite(tahoe[cell %chin% targets$cell, .(cell,drug,source_drug,dose,dose_unit,observation_id)],
  file.path(out, 'resource_Tahoe_matching_profiles.tsv'), sep = '\t')

# Exact relevant combination coverage, including a named MEKi substitution.
pair_targets <- rbindlist(list(
  CJ(cell = c('MCF7','LNCAP'), A = 'tamoxifen', B = c('mefloquine','withaferin')),
  CJ(cell = c('MCF7','LNCAP'), A = 'mefloquine', B = 'withaferin'),
  CJ(cell = c('H1792','H2009','H358'), A = 'trametinib', B = c('lestaurtinib','midostaurin')),
  CJ(cell = c('HCT116','DLD1','CACO2'), A = c('selumetinib','trametinib'), B = c('pictilisib','jnkin8'))))
pair_targets[, Jaaks_rows := vapply(seq_len(.N), function(i)
  jaaks[cell == pair_targets$cell[i] & ((drug_A == A[i] & drug_B == B[i]) | (drug_A == B[i] & drug_B == A[i])), .N], integer(1))]
pair_targets[, ALMANAC_rows := vapply(seq_len(.N), function(i)
  almanac[cell == pair_targets$cell[i] & ((drug_A == A[i] & drug_B == B[i]) | (drug_A == B[i] & drug_B == A[i])), .N], integer(1))]
fwrite(pair_targets, file.path(out, 'resource_exact_combination_overlap.tsv'), sep = '\t')
fwrite(jaaks[cell %chin% targets$cell], file.path(out, 'resource_Jaaks_matching_pairs.tsv.gz'), sep = '\t')
fwrite(almanac[cell %chin% targets$cell], file.path(out, 'resource_ALMANAC_matching_rows.tsv.gz'), sep = '\t')
# Broader transfer: the reference study uses selumetinib, whereas Jaaks tests
# trametinib + the same PI3K inhibitor. This is a named drug substitution, not
# an exact replication of the reference treatment or its colorectal context.
transfer_pair <- jaaks[(drug_A == 'trametinib' & drug_B == 'pictilisib') |
  (drug_A == 'pictilisib' & drug_B == 'trametinib')]
transfer_cells <- transfer_pair[, .(Jaaks_fit_rows = .N), by = cell]
for (nm in c('trametinib','pictilisib')) {
  transfer_cells[, (paste0('LINCS_',nm,'_samples')) := vapply(cell, function(cl)
    lincs[cell == cl & drug == nm, .N], integer(1))]
  transfer_cells[, (paste0('Tahoe_',nm,'_profiles')) := vapply(cell, function(cl)
    tahoe[cell == cl & drug == nm, .N], integer(1))]
}
transfer_cells[, reference_difference := 'Selumetinib to trametinib; breast cells rather than reference colorectal cells']
fwrite(transfer_cells, file.path(out, 'resource_MEK_PI3K_transfer_cells.tsv'), sep = '\t')
ligands <- fread('/work/data/01_CytoSig_data/meta_info.txt', skip = 1, header = FALSE,
  col.names = c('profile_id','Treatment','Condition','Duration','Dose'))
fwrite(ligands[Treatment %chin% c('HGF','EGF')], file.path(out, 'resource_HGF_EGF_profiles.tsv'), sep = '\t')
inputs <- c(tahoe_catalog_path, panel_path, prism_path, lincs_path, jaaks_paths, almanac_names_path, almanac_path)
fwrite(data.table(path = inputs, bytes = file.info(inputs)$size), file.path(out, 'resource_overlap_inputs.tsv'), sep = '\t')
writeLines(c('Overlap scope: local pinned Tahoe/PDEX, prepared39cell panel, PRISM secondary-screen, LINCS GSE70138, full local Jaaks fitted matrices and ALMANAC.',
  'Parent-compound aliases are explicit; salt/formulation differences remain visible in source names.',
  'A matching drug-cell pair is not automatically a matching dose, time, medium or genetic state.',
  'Zero in this audit means absent from the specified local release/prepared panel, not absent from all public data.',
  'Jaaks counts are fitted combination rows, not independent experiments; ALMANAC counts are source rows, not biological replicates.',
  'No engineered MYC, ligand rescue, or drug-plus-ligand measurement is inferred from ordinary single-agent profiles.'),
  file.path(root, 'RESOURCE_OVERLAP_SCOPE.txt'))
print(coverage)
print(pair_targets[Jaaks_rows > 0 | ALMANAC_rows > 0])
