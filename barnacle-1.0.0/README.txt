Using Barnacle

NOTE: Calling any command with the option "-h" will provide a full list of arguments and options
NOTE: the first two arguments to most commands are "${lib_id}" and "${input_file}", where ${lib_id} is the ID of the library being processed and ${input_file} is usually the output file of the previous stage.

SETUP:
i) Ensure that the paths and hostname and queue values in src/barnacle.cfg are correct
ii) Run setup.py to compile the portions of Barnacle that require compilation, and download and setup the default annotations
iii) If you wish to use additional gene annotation files, create appropriate "gene feature coordinates" files with src/annotation/create_gene_feature_coords.py
  USAGE: ./src/annotation/create_gene_feature_coords.py ${input_file} ${output_dir}
  e.g. ./src/annotation/create_gene_feature_coords.py annotations/UCSC_genes_ref.txt annotations/
iv) Create a configuration file for the project
  e.g.  sample_data/sample.cfg

0) Alignment processing: barnacle.pl
  USAGE: ./src/barnacle.pl -lib ${lib_id} -lib_dir ${project_dir}/${lib_id} -config ${config_path} -identify_candidates [-cluster ${cluster_submit_hostname}]
  e.g. ./src/barnacle.pl -lib SIM06 -lib_dir sample_data/SIM06 -config sample_data/sample.cfg -identify_candidates -cluster login4
  NOTE: the -cluster option is only used if you want to submit your candidate identification jobs to a cluster, in which case you should provide it with the hostname you use to submit jobs to your cluster
  NOTE: the default read length is 75 nt; if your reads are not 75 nt, then use the -read_length option to specify the correct read length
  NOTE: if your candidate identification jobs require more or less than the default 8G of memory, you can change the resource request used when submitting them to the cluster with the -cid_memory option

You can watch the status of your candidate identification (cid) jobs with: src/alignment_processing/check_status.py
  USAGE: ./src/alignment_processing/check_status.py ${cid_cluster_dir}
  e.g. ./src/alignment_processing/check_status.py sample_data/SIM06/Assembly/abyss-1.3.2/barnacle/ver_1.0.0.0/cluster_cid

If jobs fail, wait until there are no "${lib_id}-cid" jobs left on the cluster, then resubmit:
  USAGE: ./src/alignment_processing/check_status.py ${cid_cluster_dir} --resubmit --cluster-head ${cluster_submit_hostname} [--memory ${memory}]
  e.g. ./src/alignment_processing/check_status.py sample_data/SIM06/Assembly/abyss-1.3.2/barnacle/ver_1.0.0.0/cluster_cid --resubmit --cluster-head login4
  NOTE: The --memory argument is only necessary if the jobs failed because of using more than the default: 8G

1) Grouping and initialization of annotation and support:
In the Barnacle results directory there will be a "run_support.sh" script
${barnacle_dir} = ${input_dir}/barnacle/${output_ver}
WARNING: this step sometimes requires an extremely large amount of memory
You have two options for this step:
 A) run it locally
  USAGE: ./${barnacle_dir}/run_support.sh
 B) submit it to your cluster using the Barnacle submit script:
  USAGE: ./src/utils/submit.py ${job_name} ${barnacle_dir}/run_support.sh ${cluster_submit_hostname} --memory ${memory} [--email ${email_address}]
  e.g. ./src/utils/submit.py SIM06-support sample_data/SIM06/Assembly/abyss-1.3.2/barnacle/ver_1.0.0.0/run_support.sh login4 --memory 8G
  NOTE: if you use the "--email" option with your email address and your cluster is set up for it, you will receive an email if the script is killed by the cluster (e.g. if it uses more than the requested amount of memory)

NOTE: This step should automatically run steps 2-5 and submit the read-to-contig (r2c) cluster jobs, so once the r2c jobs are done, you should be able to move on to integration in step 6

If you are running jobs on several libraries, you can create a library list file and watch the status of this step:
  USAGE: ./src/grouping/check_support_status.py ${lib_list} ${project_dir} --check-cid --check-r2c
  e.g. ./src/grouping/check_support_status.py sample_data/barnacle_list.1.0.0.0.txt sample_data/ --check-cid --check-p2g --check-r2c
  NOTE: Barnacle provides you with the information needed to create the ${lib_list} file, simply concatenate the "lib_info" files that can be found in the Barnacle output directory of each library
  e.g. cat ${project_dir}/*/Assembly/${assembly_version}/barnacle/${output_ver}/lib_info > ${project_dir}/barnacle_list.${output_ver}.txt

2) Format checking: SKIP IF STEP 1 SUCCESSFUL
  USAGE: ./src/utils/check_format.py ${lib_id} ${barnacle_dir}/1_raw_candidates/${lib_id}.barnacle.data
  e.g. ./src/utils/check_format.py SIM06 sample_data/SIM06/Assembly/abyss-1.3.2/barnacle/ver_1.0.0.0/1_raw_candidates/SIM06.barnacle.data

3) Exon-boundaries: SKIP IF STEP 1 SUCCESSFUL
  USAGE: ./src/annotation/exon_bounds.py ${lib_id} ${barnacle_dir}/2_parsable_candidates/${lib_id}.barnacle.data ${gene_annotations}
  e.g. ./src/annotation/exon_bounds.py SIM06 sample_data/SIM06/Assembly/abyss-1.3.2/barnacle/ver_1.0.0.0/2_parsable_candidates/SIM06.barnacle.data annotations/ensembl65_ref.txt
  NOTE: ${gene_annotations} should be a GTF or UCSC genePredExt formatted set of gene annotations (see annotations/setup_annotations.sh for how to create properly formatted annotation sets)

4) Repeats: SKIP IF STEP 1 SUCCESSFUL
  USAGE: ./src/annotation/repeats.py ${lib_id} ${barnacle_dir}/3_with_exon_bounds/${lib_id}.barnacle.data ${repeat_files_list}
  e.g. ./src/annotation/repeats.py SIM06 sample_data/SIM06/Assembly/abyss-1.3.2/barnacle/ver_1.0.0.0/3_with_exon_bounds/SIM06.barnacle.data annotations/hg19_all_rmsk.coords,annotations/hg19_simple_repeats.coords
  NOTE: ${repeat_files_list} should be a comma-delimited list of, e.g., repeat masker and simple repeats coordinates

5) Breakpoint genes: SKIP IF STEP 1 SUCCESSFUL
  USAGE: ./src/annotation/breakpoint_genes.py ${lib_id} ${barnacle_dir}/4_with_repeats/${lib_id}.barnacle.data ${gene_feature_coords}
  e.g. ./src/annotation/breakpoint_genes.py SIM06 sample_data/SIM06/Assembly/abyss-1.3.2/barnacle/ver_1.0.0.0/4_with_repeats/SIM06.barnacle.data annotations/ensembl65_ref.exons.introns.std_chr.bed
  NOTE: ${gene_feature_coords} should be a BED-formatted file created by create_gene_feature_coords.py (which is automatically run when you run src/setup.py)

6) Read-to-contig (r2c) support:
Submission: SKIP THIS IF STEP 1 SUCCESSFUL
  USAGE ./src/support/read_to_contig/submit.py ${lib_id} ${barnacle_dir}/2_parsable_candidates/${lib_id}.barnacle.data ${input_dir}/reads_to_contigs/${lib_id}-reads-to-contigs.bam --cluster-head ${cluster_submit_hostname}
  e.g. ./src/support/read_to_contig/submit.py SIM06 sample_data/SIM06/Assembly/abyss-1.3.2/barnacle/ver_1.0.0.0/2_parsable_candidates/SIM06.barnacle.data sample_data/SIM06/Assembly/abyss-1.3.2/reads_to_contigs/SIM06-reads-to-contigs.bam --cluster-head login4
  NOTE: the default read length is 75 nt; if your reads are not 75 nt, then use the -read-length option to specify the correct read length

Watch read-to-contig status support with: support/read_to_contig/check_status.py
  e.g. ./src/support/read_to_contig/check_status.py sample_data/SIM06/Assembly/abyss-1.3.2/barnacle/ver_1.0.0.0/cluster_r2c

If jobs fail, wait until there are no "${lib_id}-r2c" jobs left on the cluster, then resubmit with:
  USAGE ./src/support/read_to_contig/check_status.py ${jobs_dir} --resubmit --cluster-head ${cluster_submit_hostname} [--memory ${memory}]
  e.g. ./src/support/read_to_contig/check_status.py sample_data/SIM06/Assembly/abyss-1.3.2/barnacle/ver_1.0.0.0/cluster_r2c --resubmit --cluster-head login4
  NOTE: The --memory argument is only necessary if the jobs failed because of using more than the default: 5G

Once all read-to-contig support jobs are complete, integrate the results:
  USAGE ./src/support/read_to_contig/integrate.py ${lib_id} ${barnacle_dir}/5_breakpoint_genes/${lib_id}.barnacle.data
  e.g. ./src/support/read_to_contig/integrate.py SIM06 sample_data/SIM06/Assembly/abyss-1.3.2/barnacle/ver_1.0.0.0/6_with_p2g/SIM06.barnacle.data

7) Filtering: src/filter/filter_groups.py
Use the "-h" option to see descriptions of the many different filtering options available.
  USAGE ./src/filter/filter_groups.py ${lib_id} ${barnacle_dir}/6_with_r2c/${lib_id}.barnacle.data
  e.g. ./src/filter/filter_groups.py SIM06 sample_data/SIM06/Assembly/abyss-1.3.2/barnacle/ver_1.0.0.0/7_with_r2c/SIM06.barnacle.data

8) Event prediction: src/prediction/predict_events.py
Again, there are many different options (use the "-h" option to see them all).
  USAGE ./src/prediction/predict_events.py ${lib_id} ${barnacle_dir}/7_filtered/${lib_id}.barnacle.pass --transcript-annotations ${gene_annotations} --transcript-sequences ${transcript_seqs} --contig-sequences ${contig_seqs}
  e.g. ./src/prediction/predict_events.py SIM06 sample_data/SIM06/Assembly/abyss-1.3.2/barnacle/ver_1.0.0.0/8_filtered/SIM06.barnacle.pass --transcript-annotations annotations/ensembl65_ref.txt --transcript-sequences annotations/Homo_sapiens.GRCh37.65.cdna.all.fa --contig-sequences sample_data/SIM06/Assembly/abyss-1.3.2/barnacle/ver_1.0.0.0/1_raw_candidates/SIM06.barnacle.contigs

9) Relative coverage: src/

DONE
