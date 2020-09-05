include_recipe 'mountetna::default'

mountetna_etna_app 'polyphemus' do
  extra_yml_config(
      metis: {
          release_bucket: node['waiver']['release_bucket'],
          restrict_bucket: node['waiver']['restrict_bucket'],
      },
      polyphemus: {
          token: node['polyphemus']['janus_token'],
      },
  )
end


