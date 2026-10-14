import React, {useState, useCallback, useEffect} from 'react';

import Checkbox from '@material-ui/core/Checkbox';
import FormControlLabel from '@material-ui/core/FormControlLabel';

import {getModels} from 'etna-js/api/magma_api';
import {getProjects} from 'etna-js/api/janus-api';

import ModelActionsModal, { ModelModalParams } from './model_actions_modal';
import {Attribute} from '../../api/magma_api';
import {SNAKE_CASE, COMMA_SEP, COMMA_SEP_WITH_SPACES, VALIDATION_TYPES} from '../../utils/edit_map';
import {ShrinkingLabelTextField} from './shrinking_label_text_field';
import ModalSelect from './modal_select';

const ontologyOf = (attribute: Attribute) => ({
  project_name: attribute.validation?.project_name || '',
  model_name: attribute.validation?.model_name || '',
  attribute_name: attribute.validation?.attribute_name || ''
});

export default function EditAttributeModal({
  onSave,
  open,
  onClose,
  attribute
}: ModelModalParams & { attribute: Attribute }) {
  const [updatedAttribute, setUpdatedAttribute] = useState({...attribute});
  const [validationType, setValidationType] = useState(
    attribute.validation ? attribute.validation.type : ''
  );
  const [validationValue, setValidationValue] = useState(
    attribute.validation ? attribute.validation.value : ''
  );
  const isArrayValidation = 'Array' === validationType;
  const isOntologyValidation = 'Ontology' === validationType;
  const [ontology, setOntology] = useState(ontologyOf(attribute));
  const [ontologyProjects, setOntologyProjects] = useState<string[]>([]);
  const [ontologyModels, setOntologyModels] = useState<{[name: string]: any}>({});

  useEffect(() => {
    if (!open || !isOntologyValidation) return;

    getProjects()
      .then((projects: {project_name: string, resource: boolean}[]) =>
        setOntologyProjects(
          projects.filter((project) => project.resource).map((project) => project.project_name)
        )
      )
      .catch(() => setOntologyProjects([]));
  }, [open, isOntologyValidation]);

  useEffect(() => {
    if (!open || !isOntologyValidation || !ontology.project_name) return;

    getModels(ontology.project_name)
      .then(({models: {project, ...models}}: {models: {[name: string]: any}}) =>
        setOntologyModels(models)
      )
      .catch(() => setOntologyModels({}));
  }, [open, isOntologyValidation, ontology.project_name]);

  const ontologyAttributes = Object.keys(
    ontologyModels[ontology.model_name]?.template.attributes || {}
  );

  const handleOnSave = useCallback(() => {
    let params = {
      ...updatedAttribute,
      attribute_name: attribute.attribute_name,
      new_attribute_name: updatedAttribute.attribute_name
    };

    if (isOntologyValidation) {
      params.validation = ontology.attribute_name ? {type: validationType, ...ontology} : null;
    } else if (validationType && validationValue) {
      params.validation = {
        type: validationType,
        value:
          isArrayValidation && !Array.isArray(validationValue)
            ? validationValue.split(',')
              // ignore leading and trailing white-space per entry
              .map((s: string) => s.trim())
            : validationValue
      };
    } else {
      params.validation = null;
    }

    onSave(params);
  }, [
    attribute,
    updatedAttribute,
    validationType,
    validationValue,
    isArrayValidation,
    isOntologyValidation,
    ontology,
    onSave
  ]);

  const reset = useCallback(() => {
    setUpdatedAttribute({...attribute});
    setValidationType( attribute.validation ? attribute.validation.type : '');
    setValidationValue( attribute.validation ? attribute.validation.value : '');
    setOntology(ontologyOf(attribute));
  }, [attribute]);

  useEffect(() => {
    reset()
  }, [attribute])

  const disabled = !updatedAttribute.attribute_name

  const handleOnCancel = useCallback(() => {
    onClose();
    reset();
  }, [onClose, reset]);

  const updateAttribute = useCallback(
    (updatePairs: [string, string | boolean][]) => {
      let tmpAttribute: Attribute = {...updatedAttribute};
      updatePairs.forEach(([key, value]) => {
        (
          tmpAttribute as Record<
            keyof Attribute,
            string | boolean | {[key: string]: any}
          >
        )[key as keyof Attribute] = value;
      });

      setUpdatedAttribute(tmpAttribute);
    },
    [updatedAttribute]
  );

  const validation_length_limit = 2000;

  return (
    <ModelActionsModal onClose={handleOnCancel} open={open} onSave={handleOnSave} title='Edit Attribute' saveDisabled={disabled}>
        <ShrinkingLabelTextField
          id='edit-attribute-name'
          label='Name (snake_case)'
          value={updatedAttribute.attribute_name}
          onChange={(e: React.ChangeEvent<any>) =>
            updateAttribute([
              ['attribute_name', e.target.value],
              ['new_attribute_name', e.target.value]
            ])
          }
          pattern={SNAKE_CASE}
        />
        <ShrinkingLabelTextField
          id='edit-attribute-description'
          label='Description'
          value={updatedAttribute.description}
          onChange={(e: React.ChangeEvent<any>) =>
            updateAttribute([['description', e.target.value]])
          }
        />
        <ShrinkingLabelTextField
          id='edit-attribute-group'
          label='Group (comma-separated list)'
          value={updatedAttribute.attribute_group}
          onChange={(e: React.ChangeEvent<any>) =>
            updateAttribute([['attribute_group', e.target.value]])
          }
          pattern={COMMA_SEP}
        />
        <ShrinkingLabelTextField
          id='edit-attribute-display-name'
          label='Display Name'
          value={updatedAttribute.display_name}
          onChange={(e: React.ChangeEvent<any>) =>
            updateAttribute([['display_name', e.target.value]])
          }
        />
        <ShrinkingLabelTextField
          id='edit-attribute-format-hint'
          label='Format Hint'
          value={updatedAttribute.format_hint}
          onChange={(e: React.ChangeEvent<any>) =>
            updateAttribute([['format_hint', e.target.value]])
          }
        />
        {attribute.attribute_type!='matrix' && (
          <ModalSelect
            id='edit-attribute-validation-type'
            label='Validation Type'
            value={validationType}
            onChange={(value: string) => {
              setValidationValue('');
              setValidationType(value);
            }}
            options={VALIDATION_TYPES}
          />
        )}
        {isOntologyValidation && (
          <>
            <ModalSelect
              id='edit-attribute-ontology-project'
              label='Ontology Project'
              value={ontology.project_name}
              onChange={(value: string) =>
                setOntology({project_name: value, model_name: '', attribute_name: ''})
              }
              options={ontologyProjects}
            />
            <ModalSelect
              id='edit-attribute-ontology-model'
              label='Ontology Table'
              value={ontology.model_name}
              onChange={(value: string) =>
                setOntology({...ontology, model_name: value, attribute_name: ''})
              }
              options={Object.keys(ontologyModels)}
            />
            <ModalSelect
              id='edit-attribute-ontology-attribute'
              label='Ontology Column'
              value={ontology.attribute_name}
              onChange={(value: string) =>
                setOntology({...ontology, attribute_name: value})
              }
              options={ontologyAttributes}
            />
          </>
        )}
        {validationType && !isOntologyValidation && (
          <ShrinkingLabelTextField
            id='edit-attribute-validation-value'
            label={`Validation ${
              isArrayValidation ?
                attribute.attribute_type=='matrix' ?
                  `Array (comma-separated list) -- UI length limit: ${validationValue.length} / 2000 characters` :
                  'Array (comma-separated list)' :
                'Regex'
            }`}
            disabled={validationValue.length > validation_length_limit}
            value={validationValue.length <= validation_length_limit ? validationValue : 'Validation too long to edit here'}
            onChange={(e: React.ChangeEvent<any>) => {
              if (e.target.value.length <= validation_length_limit) setValidationValue(e.target.value)
            }}
            pattern={isArrayValidation ? COMMA_SEP_WITH_SPACES : null}
          />
        )}
        <FormControlLabel
          control={
            <Checkbox
              onChange={(e: React.ChangeEvent<any>) =>
                updateAttribute([['hidden', e.target.checked]])
              }
              checked={updatedAttribute.hidden}
              inputProps={{'aria-label': 'controlled'}}
            />
          }
          label='Hidden'
        />
        <FormControlLabel
          control={
            <Checkbox
              onChange={(e: React.ChangeEvent<any>) =>
                updateAttribute([['read_only', e.target.checked]])
              }
              checked={updatedAttribute.read_only}
              inputProps={{'aria-label': 'controlled'}}
            />
          }
          label='Read-only'
        />
        <FormControlLabel
          control={
            <Checkbox
              onChange={(e: React.ChangeEvent<any>) =>
                updateAttribute([['restricted', e.target.checked]])
              }
              checked={updatedAttribute.restricted}
              inputProps={{'aria-label': 'controlled'}}
            />
          }
          label='Restricted'
        />
        {CONFIG.project_name === 'coprojects_template' && (
          <FormControlLabel
            control={
              <Checkbox
                onChange={(e: React.ChangeEvent<any>) =>
                  updateAttribute([['required_on_child_models', e.target.checked]])
                }
                checked={!!updatedAttribute.required_on_child_models}
                inputProps={{'aria-label': 'controlled'}}
              />
            }
            label='Required on child models'
          />
        )}
    </ModelActionsModal>
  );
}
