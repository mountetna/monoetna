class Vulcan
  # This is the temporary glue for running synchronous workflows from a given session and its inputs.
  # Most of this logic will evolve outside of vulcan depending on where we land on orchestration.
  class Orchestration
    def self.find_runnable_steps(session, workflow)
      workflow.unique_paths.map do |path|
        path.find do |step_name|
          step = workflow.find_step(step_name)
          input_files = self.input_files_for(step, workflow, session)
          output_files = self.output_files_for(step, workflow, session)

          input_files.all? { |input_file| ::File.exists?(input_file.data_path) } &&
              output_files.any? { |output_file| ! ::File.exists?(output_file.data_path) }
        end
      end
    end

    def self.input_files_for(step, workflow, session)
      step.in.map { |input_name, input_ref|
        self.create_input_storage_file(input_name, input_ref, workflow, session) }
    end

    def self.output_files_for(step, workflow, session)
      step.out.map { |output_name|
        self.create_output_storage_file(step_name, output_name, workflow, session) }
    end

    def self.cell_hash(step_name, workflow, session)
      if step_name == :primary_inputs
        raw_hash = session.primary_input[output_name]

        return raw_hash && Storage.cell_hash(
            project_name: session.project_name,
            input_files: [],
            output_filenames: [output_name],
            session_key: session.key,
            script_or_raw_hash: raw_hash,
        )
      end

      step = workflow.find_step(step_name)
      return nil if step.nil?

      input_files = []
      step.in.each do |input_logical_name, ref|
        inner_step, inner_output_name = ref
        ch = cell_hash(inner_step, workflow, session)
        return nil if ch.nil?

        input_files << Storage::StorageFile.new(
            project_name: session.project_name, cell_hash: ch,
            data_filename: inner_output_name, logical_name: input_logical_name
        )
      end

      Storage.cell_hash(
          project_name: session.project_name,
          input_files: input_files,
          output_filenames: step.out,
          session_key: session.key,
          script_or_raw_hash: step.script,
      )
    end

    def self.create_input_storage_file(input_name, reference, workflow, session)
      step_name, output_file = reference
      ch = cell_hash(step_name, workflow, session)

      return nil if ch.nil?

      Storage::StorageFile.new(
          project_name: session.project_name, cell_hash: ch,
          data_filename: output_file, logical_name: input_name
      )
    end

    def self.create_output_storage_file(step_name, output_name, workflow, session)
      ch = cell_hash(step_name, workflow, session)

      return nil if ch.nil?

      Storage::StorageFile.new(
          project_name: session.project_name, cell_hash: ch,
          data_filename: output_name, logical_name: output_name,
      )
    end
  end
end