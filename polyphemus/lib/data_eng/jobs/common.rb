def build_pipeline_state_dir(path_to_write_files, run_id)
  dir_path = File.join(path_to_write_files, run_id)
  FileUtils.mkdir_p(dir_path)
  dir_path
end

def remove_dscolab_prefix(path)
    path.gsub('DSCOLAB_', '')
end