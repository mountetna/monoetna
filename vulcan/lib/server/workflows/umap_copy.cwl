cwlVersion: v1.1
class: Workflow

inputs:
  6_Cluster_Differential_Expression__dge_method:
    type: string
    default: 'wilcoxon'
    label: 'testing method'
    doc: 'A string indicating what scanpy method option to use for calculating differential expression. Options are: "logreg", "t-test", "wilcoxon", "t-test_overestim_var". See documentation for "scanpy.tl.rank_genes_groups" for further details.'

outputs:
  the_data:
    type: File
    outputSource: thing/hash

steps:
  thing:
    run: scripts/make_hash.cwl
    label: 'thing'
    in: 
      a: 6_Cluster_Differential_Expression__dge_method
    out: [hash]
  makeThing:
    run: ui-queries/multiple-string.cwl
    label: 'do the thing'
    in:
      a: thing/hash
    out: [out]
