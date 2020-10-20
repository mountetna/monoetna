require 'etna'

class IpiAddFlowModelMigration

  def initialize(magma_client:, magma_crud:)
    @magma_client = magma_client
    @magma_crud = magma_crud
    @project_name = 'ipi'

    @flow_stains = ['dc', 'innate', 'nktb', 'sort', 'treg']
    @flow_attribute_suffixes = ['file', 'flag', 'notes']

    @suffix_to_new_attr_map = {
      file: 'fcs_file',
      flag: 'quality_flag',
      notes: 'notes'
    }

    @samples = fetch_sample_documents
  end

  def sample_flow_attributes
    @flow_stains.map { |stain|
      @flow_attribute_suffixes.map { |suffix|
        "#{stain}_#{suffix}"
      }
    }.flatten
  end

  def create_flow_model
    create_flow_request = Etna::Clients::Magma::UpdateModelRequest.new(
      project_name: @project_name,
      actions: [Etna::Clients::Magma::AddModelAction.new(
        model_name: 'flow',
        parent_model_name: 'sample',
        parent_link_type: 'collection',
        identifier: 'stain_name'
      )]
    )
    @magma_client.update_model(create_flow_request)

    add_attributes_request = Etna::Clients::Magma::UpdateModelRequest.new(
      project_name: @project_name,
      actions: [Etna::Clients::Magma::AddAttributeAction.new(
        model_name: 'flow',
        type: 'string',
        attribute_name: 'stain',
        description: 'Flow Cytometry stain panel.',
        display_name: 'Stain',
        validation: {
          "type": "Array",
          "value": ["dc", "innate", "nktb", "sort", "treg"]
        }
      ), Etna::Clients::Magma::AddAttributeAction.new(
        model_name: 'flow',
        type: 'file',
        attribute_name: 'fcs_file',
        description: 'FACS format file.',
        display_name: 'FCS file'
      ), Etna::Clients::Magma::AddAttributeAction.new(
        model_name: 'flow',
        type: 'string',
        attribute_name: 'quality_flag',
        description: 'Quality flag for the stain.',
        display_name: 'Qualify flag'
      ), Etna::Clients::Magma::AddAttributeAction.new(
        model_name: 'flow',
        type: 'string',
        attribute_name: 'notes',
        description: 'Notes for issues with the stain.',
        display_name: 'Notes'
      )]
    )

    @magma_client.update_model(add_attributes_request)
  end

  def create_population_model
    create_population_request = Etna::Clients::Magma::UpdateModelRequest.new(
      project_name: @project_name,
      actions: [Etna::Clients::Magma::AddModelAction.new(
        model_name: 'population',
        parent_model_name: 'flow',
        parent_link_type: 'table'
      )]
    )
    @magma_client.update_model(create_population_request)

    add_attributes_request = Etna::Clients::Magma::UpdateModelRequest.new(
      project_name: @project_name,
      actions: [Etna::Clients::Magma::AddAttributeAction.new(
        model_name: 'population',
        type: 'string',
        attribute_name: 'name',
        description: 'Name of this population.'
      ), Etna::Clients::Magma::AddAttributeAction.new(
        model_name: 'population',
        type: 'integer',
        attribute_name: 'count',
        description: 'Number of cells.'
      ), Etna::Clients::Magma::AddAttributeAction.new(
        model_name: 'population',
        type: 'string',
        attribute_name: 'ancestry',
        description: 'Chain of parent populations.'
      )]
    )

    @magma_client.update_model(add_attributes_request)
  end

  def create_mfi_model
    create_mfi_request = Etna::Clients::Magma::UpdateModelRequest.new(
      project_name: @project_name,
      actions: [Etna::Clients::Magma::AddModelAction.new(
        model_name: 'mfi',
        parent_model_name: 'population',
        parent_link_type: 'table'
      )]
    )
    @magma_client.update_model(create_mfi_request)

    add_attributes_request = Etna::Clients::Magma::UpdateModelRequest.new(
      project_name: @project_name,
      actions: [Etna::Clients::Magma::AddAttributeAction.new(
        model_name: 'mfi',
        type: 'string',
        attribute_name: 'name',
        description: 'Name of antibody.'
      ), Etna::Clients::Magma::AddAttributeAction.new(
        model_name: 'mfi',
        type: 'string',
        attribute_name: 'fluor',
        description: 'Color of this channel.'
      ), Etna::Clients::Magma::AddAttributeAction.new(
        model_name: 'mfi',
        type: 'float',
        attribute_name: 'value',
        description: 'Mean fluorescence intensity.'
      )]
    )

    @magma_client.update_model(add_attributes_request)
  end

  def copy_flow_data
    # Copy the existing file, quality_flag, and notes for each stain
    #   from the Sample model to the new Flow model

    # Do the updates per Sample, so we don't exceed the Rack key space
    #   on update
    @samples.document_keys.each { |sample_name|
      sample = @samples.document(sample_name)

      update_records_request = Etna::Clients::Magma::UpdateRequest.new(project_name: @project_name)

      @flow_stains.each { |stain|
        @flow_attribute_suffixes.each { |suffix|
          old_attribute_name = "#{stain}_#{suffix}"
          new_attribute_name = @suffix_to_new_attr_map[suffix.to_sym]
          flow_record_name = "#{sample['sample_name']}.flow.#{stain}"

          # Don't bother copying nil values to the Flow model
          next if sample[old_attribute_name] == nil

          update_records_request.send(
            :update_revision,
            'flow',
            flow_record_name,
            sample: sample_name,
            stain: stain,
            "#{new_attribute_name}": format_attribute_value(suffix, sample[old_attribute_name]))
        }
      }

      # Only send the update request if there are revisions
      @magma_client.update(update_records_request) if update_records_request.revisions.length > 0
    }
  end

  def hide_sample_columns
    sample_flow_attributes.each { |attribute_name|
      hide_columns_request = Etna::Clients::Magma::UpdateModelRequest.new(project_name: @project_name)

      hide_columns_request.add_action(
        Etna::Clients::Magma::UpdateAttributeAction.new(
          model_name: 'sample',
          attribute_name: attribute_name,
          hidden: true
        )
      )

      @magma_client.update_model(hide_columns_request)
    }
  end

  def format_attribute_value(suffix, value)
    return value unless suffix == 'file'

    # So this is kind of wonky ... taking value['path']
    #   as the name of the file on Metis works because
    #   the source file is in the Magma bucket, which
    #   was named according to the value['path'] file name.
    # So even though we don't know the original Metis
    #   path or file name, we can still create a copy of the
    #   Magma copy of the file.
    return {
      path: "metis://#{@project_name}/magma/#{value['path']}"
    } if value&.dig('path')

    return {
      path: nil
    }
  end

  def fetch_sample_documents
    request = Etna::Clients::Magma::RetrievalRequest.new(project_name: @project_name)
    request.model_name = 'sample'
    request.attribute_names = sample_flow_attributes
    request.record_names = 'all'
    @magma_client.retrieve(request).models.model('sample').documents
  end

  def execute
    create_flow_model
    copy_flow_data
    hide_sample_columns
    create_population_model
    create_mfi_model
  end
end