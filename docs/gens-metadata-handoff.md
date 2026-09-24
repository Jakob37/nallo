# Gens metadata implementation handoff

## Decisions from follow-up

- Generate coverage and ROH metadata for every sample when Gens metadata is
  enabled. Use a dedicated metadata flag so existing Gens coverage/BAF input
  generation remains independently configurable.
- Generate UPD metadata and tracks only for children whose UPD analysis actually
  ran. Parents, singletons, and runs with UPD disabled should retain their
  coverage and ROH metadata without UPD values. A completed UPD analysis with
  no sites or regions is distinct from an analysis that did not run.
- Follow the current OLWGS metadata script's `Non-informative` calculation for
  now: sum `UPD_PATERNAL_ORIGIN`, `UPD_MATERNAL_ORIGIN`, and `ANTI_UPD`.
  **TODO:** Review this label and calculation later. The UPD producer calls
  these informative categories and also emits `UNINFORMATIVE`; the older OLWGS
  `bin/upd_table.pl` uses a different calculation.

## Original investigation

The user wants Gens metadata calculations in Nallo similar to the adjacent
`/home/jakob/src/olwgs` (nextflow_wgs) pipeline. The initial investigation
covered both repositories. The implementation on `jw/gens-metadata` now adds a
metadata formatter and workflow wiring; full pipeline validation remains open.

Nallo workspace: `/home/jakob/src/nallo`. OLWGS is readable at
`/home/jakob/src/olwgs`. The working tree was clean before this note was added.
Read any applicable repository instructions before editing.

## Reference implementation in OLWGS

- `bin/prepare_gens_v4_input.py`: Python standard-library implementation of
  metadata calculations and track conversion.
- `main.nf:1876`: `generate_gens_v4_meta` process and its input/output contract.
- `main.nf:503`: orchestration and assembly of metadata inputs.
- `main.nf:1710`: ROH calling, using
  `bcftools roh --rec-rate 1e-9 --AF-tag GNOMADAF`.
- `bin/create_gens_case_yaml.py`: separate packaging of files for `gens load case`.
- `tests/test_create_gens_case_yaml.py`: case YAML tests, not calculation tests.

The metadata script consumes bcftools ROH output, UPD regions and sites,
GATK copy ratios, a reference sequence dictionary, sample ID, and sex (`M`/`F`).
It produces:

| File | Contents |
| --- | --- |
| `<sample>.meta.tsv` | Header `type\tvalue`; `%Autosomal LOH` |
| `<sample>.chrom_meta.tsv` | Header `Chromosome\ttype\tvalue`; estimated copy number and UPD metrics |
| `<sample>.gens_track.roh.bed` | ROH intervals labelled LOH, with colour in the ninth field |
| `<sample>.gens_track.upd.bed` | UPD intervals labelled by parental origin, with colour in the ninth field |

Calculation details in the existing script:

- Filter ROH rows by sample and quality >=85 (configurable in script).
- Autosomal LOH percentage: sum of retained autosomal ROH lengths divided by
  summed autosomal reference lengths, multiplied by 100.
- Estimated chromosome copy number:
  `baseline_ploidy * 2 ** mean(log2_copy_ratio)`, rounded to two decimals.
  Baseline is 1 for male X/Y and 2 otherwise; female Y is omitted.
  The mean is unweighted across bins.
- OLWGS passes **denoised** ratios to metadata, although its BAF/coverage
  preparation uses standardized ratios.
- UPD statistics include total sites, mismatch mother/father, anti-UPD, and a
  field labelled non-informative, with percentages rounded to one decimal.
- `analysis_mode=single` omits UPD site statistics. OLWGS still supplies empty
  UPD files through its workflow.
- The Nextflow process only calculates metadata for `type == "proband"`;
  other samples get touched empty files.
- Case YAML attaches `meta_files` and `sample_annotations` for probands;
  UPD annotations are attached only for trio cases.

## Existing Nallo integration points

### Coverage and BAF

`subworkflows/local/prepare_gens_inputs/main.nf` already:

1. Computes coverage bins with mosdepth.
2. Formats coverage for GATK.
3. Selects a sex-specific panel of normals (`meta.sex == 1` male, `2` female).
4. Calls `GATK4_DENOISEREADCOUNTS`.
5. Sends standardized ratios and gVCFs to `PREPARECOVANDBAF`.

Both standardized and denoised outputs exist in the GATK module, but the
subworkflow currently emits only indexed BAF and coverage BED files. Expose
denoised ratios for metadata; no new coverage calculation is needed.

Relevant files:

- `modules/nf-core/gatk4/denoisereadcounts/main.nf`
- `conf/modules/prepare_gens_input.config`
- `workflows/nallo.nf:623`: current Gens invocation, before SNV annotation.
- `workflows/nallo.nf:1282`: Gens output emits.
- `main.nf`: wrapper emits, `ch_gens` assembly around line 609, workflow output
  publishing around line 790 to `gens/${meta.id}/`.

### ROH

`subworkflows/local/chromograph/main.nf` currently owns `BCFTOOLS_ROH` inside
the autozygosity branch. It passes the raw ROH file to `RHOCALL_VIZ`, but emits
only Chromograph plots.

- Existing reusable module: `modules/nf-core/bcftools/roh/main.nf`.
- Current caller settings: `conf/modules/chromograph.config`, including
  `--AF-tag ${params.bcftools_roh_af_tag}` and `--skip-indels`.
- `workflows/nallo.nf` around line 880 prepares unfiltered Echtvar-annotated
  SNVs for Chromograph/Peddy and splits to sample VCFs. Investigate reuse and
  update gating so this preparation also runs when Gens ROH metadata needs it.
- Chromograph's sample VCF metadata is reduced to `[id: sample_id]`.

### UPD

`subworkflows/local/call_upd/main.nf` already emits:

- `regions = UPD_REGIONS.out.bed`
- `sites = UPD_SITES.out.bed`

It selects every child with both named parents present in the same family,
uses unfiltered Echtvar-annotated family variants, and emits child metadata
containing `id`, `family_id`, `mother`, and `father`.

Invocation is in `workflows/nallo.nf:843`, within SNV annotation. It respects
`params.skip_upd`; it does not manufacture results for ineligible samples.

## Recommended design

1. Extract ROH calling into an independent `CALL_ROH` subworkflow. Calculate
   ROH once and share it with Chromograph and Gens. Chromograph should consume
   ROH rather than owning the calculation. Gens metadata must work with
   Chromograph disabled.
2. Expose denoised copy ratios from `PREPARE_GENS_INPUTS`.
3. Add a small local metadata module with a Python adapter based on the OLWGS
   script, plus a `PREPARE_GENS_METADATA` orchestration subworkflow if useful.
   Use normal Nallo module conventions for environment/container, versions,
   stubs, and tests. Avoid modifying vendored nf-core modules unnecessarily.
4. Assemble inputs after the required analysis outputs are available. The
   existing early Gens coverage block can stay where it is; metadata can run
   later after ROH/UPD outputs are available.
5. Join using explicit sample identity, retaining family identity where
   available. Do not join whole metadata maps: coverage, ROH, and UPD carry
   different fields. Check sample-ID uniqueness rules before choosing keys.
   Do not copy OLWGS's unkeyed `combine` chains: these can mix families or
   multiply records in multi-family runs.
6. Treat absent analyses as absent inputs. Missing UPD must not suppress
   coverage/ROH metadata for parents, singletons, or skipped-UPD runs.
   Distinguish an analysis that ran and found zero events from one not run.
   Avoid fake placeholder analysis files and misleading zero statistics.
7. Extend the workflow emits and existing `gens/{sample}/` publishing for
   metadata and tracks; update documentation and dependency validation.
8. Keep Gens case YAML/path packaging separate from numerical calculations.
   It is an optional additional integration step, not an established part of
   this implementation's scope. Do not copy OLWGS cron/deployment behavior.

Suggested policy is coverage/ROH metadata for every available sample and UPD
metadata for eligible children, leaving proband-only presentation downstream.
This is a recommendation, not a user-confirmed requirement.

Review parameter gating deliberately: current Gens preparation requires
mapping/SNV calling and reference resources, while ROH also requires suitable
AF annotations. Do not accidentally make existing coverage/BAF-only usage
require all metadata analyses. Relevant validation is in
`subworkflows/local/utils_nfcore_nallo_pipeline/main.nf`; schema/defaults/docs
must agree. Decide whether metadata gets its own opt-in/skip setting.

## Issues to resolve before porting calculations verbatim

- **Chromosome names:** OLWGS hardcodes bare `1`–`22`, `X`, `Y`. A `chr`-prefixed
  dictionary yields an empty autosomal denominator. Normalize canonical
  chromosome identity consistently across inputs while preserving whatever
  output names the Gens consumer requires. Handle other contigs deliberately.
- **Missing coverage:** the OLWGS script substitutes log2 ratio zero for a
  missing chromosome, falsely reporting normal copy number. Omit or explicitly
  represent unavailable values instead. Validate nonfinite/invalid values.
- **UPD labels:** `Non-informative` is currently the sum of
  `UPD_PATERNAL_ORIGIN`, `UPD_MATERNAL_ORIGIN`, and `ANTI_UPD`. Verify producer
  category semantics before preserving or correcting this calculation/label.
- **Coordinates:** OLWGS directly copies ROH starts to BED and computes lengths
  as `end - start`. Verify bcftools ROH and UPD coordinate conventions, convert
  to BED correctly, and test boundaries. Do not assume both producers use the
  same conventions. Consider overlapping ROH intervals when calculating totals.
- **ROH parser:** it checks for fewer than seven fields but subsequently reads
  an eighth field. Parse explicit region records and validate their fields.
- **Caller parity:** Nallo's ROH settings differ from OLWGS (`--skip-indels`
  versus explicit `--rec-rate 1e-9`). Choose and document settings rather than
  promising identical metadata merely because the formatter is reused.
- **Coverage aggregation:** the reference calculates the exponential of the
  unweighted mean log2 ratio, not an arithmetic mean copy number or a
  length-weighted estimate. Preserve intentionally for parity, or explicitly
  document and test a changed definition.
- **Reference denominator:** ensure usable autosomal lengths and fail clearly
  on missing/zero denominator. A reference FAI may be a convenient Nallo input;
  OLWGS currently reads a sequence dictionary.
- **Sex:** convert Nallo's numeric sex explicitly; do not silently treat an
  unknown value as female. Existing Gens validation requires known sex.

These are findings from source inspection; producer/consumer semantics and
long-read numerical behavior have not yet been independently validated.

## Validation plan

Use small meaningful Python fixtures to establish expected numerical values:

- Known autosomal ROH percentage, quality filtering, sample filtering, and
  coordinate conversion.
- Both `chr`-prefixed and bare chromosome inputs; missing reference lengths.
- Known log2 ratios and sex scaling; missing chromosomes and invalid values.
- Every UPD category, zero sites, and absent UPD versus an empty valid result.
- Exact TSV headers and annotation field layout expected by Gens.
- A compatibility fixture comparing intentional unchanged calculations with
  OLWGS, documenting any deliberate differences.

Use nf-test at module/subworkflow level for actual channel behavior:

- Singleton and trio; more than one family; multiple eligible siblings.
- Gens metadata with Chromograph disabled.
- Chromograph and Gens enabled together without duplicate ROH calling.
- UPD skipped or no eligible trio, without losing other metadata outputs.
- Existing coverage/BAF-only behavior and output publishing.

Update affected snapshots and run appropriate existing tests after refactoring.
Do not claim end-to-end or numerical parity until tested.
