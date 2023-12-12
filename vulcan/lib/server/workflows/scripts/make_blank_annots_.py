from archimedes.functions.dataflow import input_json, input_var,output_json

cluster_meta=input_var('cluster_meta')
discrete_metadata_summary=input_json('discrete_metadata_summary')

clusts = discrete_metadata_summary[cluster_meta]

output_json(
    dict(
        [
            str(clust),
            str(clust)
        ] for clust in range(max(
            [int(str) for str in clusts]
            )+1)
    ),
    'blank_annots.json'
)