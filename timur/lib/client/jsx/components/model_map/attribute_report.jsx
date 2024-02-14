import React, {useState, useEffect, useCallback, useMemo} from 'react';
import {useDispatch} from 'react-redux';
import Grid from '@material-ui/core/Grid';
import Card from '@material-ui/core/Card';
import CardContent from '@material-ui/core/CardContent';
import {makeStyles} from '@material-ui/core/styles';
import Button from '@material-ui/core/Button';
import IconButton from '@material-ui/core/IconButton';
import Tooltip from '@material-ui/core/Tooltip';
import CloseIcon from '@material-ui/icons/Close';
import EditIcon from '@material-ui/icons/Edit';
import DeleteIcon from '@material-ui/icons/Delete';

import {requestAnswer} from 'etna-js/actions/magma_actions';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {selectModels} from 'etna-js/selectors/magma';

import {
  updateAttribute,
  removeAttribute,
  removeLink
} from '../../api/magma_api';
import MapHeading from './map_heading';
import EditAttributeModal from './edit_attribute_modal';
import {
  EDITABLE_ATTRIBUTE_TYPES,
  REMOVABLE_ATTRIBUTE_TYPES,
  UNREMOVABLE_ATTRIBUTE_NAMES,
  UNEDITABLE_ATTRIBUTE_NAMES
} from '../../utils/edit_map';
import useMagmaActions from './use_magma_actions';

const ATT_ATTS = [
  'attribute_type',
  'link_model_name',
  'description',
  'display_name',
  'format_hint',
  'attribute_group',

  'validation_type',
  'validation',

  'restricted',
  'read_only',
  'hidden'
];

const useStyles = makeStyles((theme) => ({
  clauseTitle: {
    fontSize: '1.2rem'
  },
  attribute_report: {
    flex: '1 1 auto',
    borderTop: '1px solid #bbb'
  },
  attribute_card: {
    background: '#eee',
    padding: '15px',
    height: 'calc(100% - 30px)'
  },
  content: {
    height: 'calc(100% - 64px)',
    overflowY: 'auto'
  },
  type: {
    color: 'gray'
  },
  value: {
    color: '#131',
    borderBottom: '1px solid rgba(34, 139, 34, 0.1)',
    maxHeight: '90px',
    wordWrap: 'break-word',
    width: '0px',
    overflowY: 'auto'
  },
  button: {
    margin: '0.5rem'
  }
}));

const ManageAttributeActions = ({
  handleEditAttribute,
  attribute,
  handleRemoveAttribute,
  handleRemoveLink
}) => {
  const classes = useStyles();

  const models = useReduxState((state) => selectModels(state));

  const [ showEditAttributeModal, setShowEditAttributeModal ] = useState(false);

  const reciprocalAttribute = useMemo(() => {
    if (!attribute.link_model_name || !attribute.link_attribute_name)
      return null;

    return models[attribute.link_model_name].template.attributes[
      attribute.link_attribute_name
    ];
  }, [models, attribute]);

  const handleConfirmRemove = useCallback(() => {
    if (
      confirm(
        'Removing an attribute is not reversible -- are you sure? If you may want to use the attribute later, you can make it "hidden" so users cannot see it.'
      )
    ) {
      handleRemoveAttribute();
    }
  }, [handleRemoveAttribute]);

  const handleConfirmRemoveLink = useCallback(() => {
    if (
      confirm(
        'Removing a link is not reversible, and affects both models -- are you sure?'
      )
    ) {
      handleRemoveLink();
    }
  }, [handleRemoveLink]);

  const isEditableAttribute = useMemo(() => {
    return (
      EDITABLE_ATTRIBUTE_TYPES.includes(attribute.attribute_type) &&
      !UNEDITABLE_ATTRIBUTE_NAMES.includes(attribute.attribute_name)
    );
  }, [attribute, EDITABLE_ATTRIBUTE_TYPES, UNEDITABLE_ATTRIBUTE_NAMES]);

  const isRemovableAttribute = useMemo(() => {
    return (
      REMOVABLE_ATTRIBUTE_TYPES.includes(attribute.attribute_type) &&
      !UNREMOVABLE_ATTRIBUTE_NAMES.includes(attribute.attribute_name)
    );
  }, [attribute, REMOVABLE_ATTRIBUTE_TYPES, UNREMOVABLE_ATTRIBUTE_NAMES]);

  const isLinkAttribute = useMemo(() => {
    if (!reciprocalAttribute) return false;

    return (
      'link' === attribute.attribute_type ||
      'link' === reciprocalAttribute.attribute_type
    );
  }, [attribute, reciprocalAttribute]);

  if (isLinkAttribute)
    return (
      <Tooltip title='Remove Link'>
        <Button
          startIcon={<DeleteIcon />}
          onClick={handleConfirmRemoveLink}
          size='small'
          color='primary'
          className={classes.button}
        >
          Remove Link
        </Button>
      </Tooltip>
    );

  if (!isEditableAttribute) return null;

  return (
    <>
      {isRemovableAttribute && (
        <Tooltip title='Remove Attribute'>
          <Button
            startIcon={<DeleteIcon />}
            onClick={handleConfirmRemove}
            size='small'
            color='primary'
            className={classes.button}>Remove</Button>
        </Tooltip>
      )}
      <Tooltip title='Edit Attribute'>
        <Button
          startIcon={<EditIcon />}
          onClick={() => setShowEditAttributeModal(true)}
          size='small'
          color='secondary'
          className={classes.button}>Edit</Button>
      </Tooltip>
      <EditAttributeModal
        open={showEditAttributeModal}
        onClose={ () => setShowEditAttributeModal(false)}
        attribute={attribute}
        onSave={handleEditAttribute}/>
    </>
  );
};

const AttributeReport = ({attribute, model_name, isAdminUser, dismiss}) => {
  const dispatch = useDispatch();
  const [sample, setSample] = useState(null);

  const {executeAction} = useMagmaActions();

  const showSample = () => {
    requestAnswer({
      query: [model_name, '::distinct', attribute.attribute_name]
    })(dispatch).then(({answer}) => setSample(answer));
  };

  useEffect(() => {
    setSample(null);
  }, [attribute]);

  const classes = useStyles();

  const handleEditAttribute = useCallback(
    (params) => {
      executeAction(
        updateAttribute({
          model_name,
          ...params
        })
      );
    },
    [executeAction, model_name, updateAttribute]
  );

  const handleRemoveAttribute = useCallback(() => {
    executeAction(
      removeAttribute({
        model_name,
        attribute_name: attribute.name
      })
    );
  }, [executeAction, model_name, removeAttribute, attribute]);

  const handleRemoveLink = useCallback(() => {
    executeAction(
      removeLink({
        model_name,
        attribute_name: attribute.name
      })
    );
  }, [executeAction, model_name, removeLink, attribute]);

  const showSampleButton = useMemo(() => {
    if (!attribute?.attribute_type) return false;

    const allowedTypes = ['string', 'identifier'];
    return allowedTypes.includes(attribute.attribute_type);
  }, [attribute]);

  return (
    <Grid className={classes.attribute_report}>
      <Card className={classes.attribute_card}>
        <MapHeading name='Attribute' title={attribute.attribute_name}>
          {isAdminUser && (
            <ManageAttributeActions
              attribute={attribute}
              handleEditAttribute={handleEditAttribute}
              handleRemoveAttribute={handleRemoveAttribute}
              handleRemoveLink={handleRemoveLink}
            />
          )}
          {showSampleButton && (
            <Tooltip title='Show data sample'>
              <Button
                onClick={showSample}
                size='small'
                color='secondary'
                className={classes.button}
              >
                Sample
              </Button>
            </Tooltip>
          )}
          <IconButton size='small' onClick={dismiss}><CloseIcon/></IconButton>
        </MapHeading>
        <CardContent className={classes.content}>
          {ATT_ATTS.map((att) => {
            switch (att) {
              case 'validation':
                return [
                  att,
                  attribute.validation
                    ? JSON.stringify(attribute.validation.value)
                    : null
                ];
              case 'validation_type':
                return [att, attribute.validation?.type];
              default:
                return [att, attribute[att]];
            }
          })
            .filter(([name, value]) => value)
            .map(([name, value]) => (
              <Grid container key={name}>
                <Grid item xs={3} className={classes.type}>
                  {name}
                </Grid>
                <Grid item xs={9} className={classes.value}>
                  {value}
                </Grid>
              </Grid>
            ))}
          {sample && (
            <Grid container>
              <Grid item xs={3} className={classes.type}>
                sample
              </Grid>
              <Grid item xs={9} className={classes.value}>
                {sample.length > 0 ? JSON.stringify(sample) : <i>No values</i>}
              </Grid>
            </Grid>
          )}
        </CardContent>
      </Card>
    </Grid>
  );
};

export default AttributeReport;
