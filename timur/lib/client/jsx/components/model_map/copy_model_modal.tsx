import React, {useState, useCallback, useEffect, useMemo} from 'react';

import {makeStyles} from '@material-ui/core/styles';
import Typography from '@material-ui/core/Typography';
import Grid from '@material-ui/core/Grid';

import {selectModels} from 'etna-js/selectors/magma';
import {useReduxState} from 'etna-js/hooks/useReduxState';

import {EDITABLE_OPTIONS} from '../../utils/edit_map';
import {SelectProjectModel} from '../select_project_model';
import ModelActionsModal, { ModelModalParams } from './model_actions_modal';

import {getDocuments} from 'etna-js/api/magma_api';
import {Template} from '../../api/magma_api';
import ModelAttributesTable from './model_attributes_table';
import {isLink} from '../../utils/attributes';

const useStyles = makeStyles((theme) => ({
  container: {
    overflow: 'hidden',
    height: '100%',
    width: '100%'
  },
  switch: {
    marginTop: '1rem'
  },
  filler: {
    height: '48px'
  },
  attributes: {
    width: '100%'
  },
  selected_attribute_text: {
    margin: '10px 0px',
    padding: '5px 14px',
    border: '1px solid #ccc',
    borderRadius: '2px',
    background: '#eee'
  },
  attribute_text: {
    margin: '10px 0px',
    padding: '5px 14px',
    border: '1px solid #ccc',
    borderRadius: '2px'
  }
}));

export default function CopyModelModal({
  modelName,
  open,
  onSave,
  onClose
}: ModelModalParams & { modelName: string }) {
  const [templateProjectName, setTemplateProjectName] = useState(null);
  const [templateModelName, setTemplateModelName] = useState('');
  const [template, setTemplate] = useState<Template|null>(null);
  const [error, setError] = useState(null);
  const [selected, setSelected] = useState({});

  const classes = useStyles();

  const models = useReduxState((state: any) => selectModels(state));

  const currentModel = models[modelName];

  useEffect( () => {
    if (templateModelName && templateProjectName) {
      getDocuments(
        {
          project_name: templateProjectName,
          model_name: templateModelName,
          record_names: [],
          attribute_names: 'all'
        },
        fetch
      ).then(({models}) => {
        const newTemplate = models[templateModelName].template;

        Object.keys(newTemplate.attributes).filter(
          att_name => {
            const attribute = newTemplate.attributes[att_name];
            return isLink(attribute)
              || attribute.hidden
              || attribute.attribute_type == 'identifier'
              || (att_name in currentModel.template.attributes)
          }
        ).forEach(
          att_name => delete newTemplate.attributes[att_name]
        );

        setTemplate(newTemplate);
      });
    }
  }, [ templateModelName ]);

  const setTemplateModel = useCallback( newModelName => {
    setTemplateModelName(newModelName);
    setTemplate(null);
    setError(null);
    setSelected({});
  }, [] );

  const setTemplateProject = useCallback( newProjectName => {
    setTemplateProjectName(newProjectName);
    setTemplateModel('');
  }, []);

  const handleOnSave = useCallback(() => {
    if (!template) return;

    const attributes = Object.keys(selected).map(
      attribute_name => {
        const attribute = template.attributes[attribute_name];

        const filtered_attribute = Object.fromEntries(
          Object.entries(attribute)
            .filter(([key]) => EDITABLE_OPTIONS.includes(key))
        );

        return {
          ...filtered_attribute,
          type: attribute.attribute_type,
          attribute_name,
          model_name: modelName
        };
      }
    );
    onSave(attributes);
  }, [selected, template]);

  const reset = useCallback(() => {
    setTemplateProjectName(null);
    setTemplateModelName('');
    setTemplate(null);
    setError(null);
    setSelected({});
  }, [])

  const handleOnCancel = useCallback(() => {
    onClose();
    reset();
  }, []);

  const numSelected = Object.keys(selected).length;

  const numAttributes = template?.attributes ? Object.keys(template.attributes).length : 0;

  const disabled = numSelected == 0;

  return (
    <ModelActionsModal onClose={handleOnCancel} open={open} onSave={handleOnSave} title='Copy Model' saveDisabled={disabled}>
      <Grid container direction='column' className={classes.container}>
      <SelectProjectModel
        project_name={templateProjectName}
        setProjectName={setTemplateProject}
        model_name={templateModelName}
        setModelName={setTemplateModel}
        error={error}
        setError={setError}
      />
      { template && <>
        { numSelected > 0
            ? <Typography className={classes.selected_attribute_text}>{numSelected} selected</Typography>
            : <Typography className={classes.attribute_text}>{numAttributes > 0 ? 'Select' : 'No'} attributes to copy (links, identifiers and overlapping names are filtered automatically)</Typography>
        }
        { numAttributes > 0 && <ModelAttributesTable
            className={ classes.attributes }
            template={template}
            selected={selected}
            setSelected={setSelected}
          />
        }
      </>}
    </Grid>
      
    </ModelActionsModal>
  );
}
