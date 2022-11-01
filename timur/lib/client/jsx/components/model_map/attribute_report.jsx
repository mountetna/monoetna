import React, {useState, useEffect, useCallback, useMemo} from 'react';
import {useDispatch} from 'react-redux';
import Grid from '@material-ui/core/Grid';
import Card from '@material-ui/core/Card';
import CardContent from '@material-ui/core/CardContent';
import {makeStyles} from '@material-ui/core/styles';
import Button from '@material-ui/core/Button';
import Tooltip from '@material-ui/core/Tooltip';
import EditIcon from '@material-ui/icons/Edit';
import DeleteIcon from '@material-ui/icons/Delete';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {requestAnswer} from 'etna-js/actions/magma_actions';
import {useModal} from 'etna-js/components/ModalDialogContainer';
import {selectUser} from 'etna-js/selectors/user-selector';
import {isAdmin} from 'etna-js/utils/janus';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {selectModels} from 'etna-js/selectors/magma';

import {showMessages} from 'etna-js/actions/message_actions';
import {addTemplatesAndDocuments} from 'etna-js/actions/magma_actions';
import {
  updateAttribute,
  removeAttribute,
  removeLink
} from '../../api/magma_api';
import MapHeading from './map_heading';
import EditAttributeModal from './edit_attribute_modal';
import {EDITABLE_ATTRIBUTE_TYPES} from '../../utils/edit_map';

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
    height: '50%',
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
  const {openModal} = useModal();

  const models = useReduxState((state) => selectModels(state));

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
    return EDITABLE_ATTRIBUTE_TYPES.includes(attribute.attribute_type);
  }, [attribute, EDITABLE_ATTRIBUTE_TYPES]);

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
      <Tooltip title='Remove Attribute'>
        <Button
          startIcon={<DeleteIcon />}
          onClick={handleConfirmRemove}
          size='small'
          color='primary'
          className={classes.button}
        >
          Remove
        </Button>
      </Tooltip>
      <Tooltip title='Edit Attribute'>
        <Button
          startIcon={<EditIcon />}
          onClick={() => {
            openModal(
              <EditAttributeModal
                attribute={attribute}
                onSave={handleEditAttribute}
              />
            );
          }}
          size='small'
          color='secondary'
          className={classes.button}
        >
          Edit
        </Button>
      </Tooltip>
    </>
  );
};

const AttributeReport = ({attribute, model_name, counts}) => {
  const dispatch = useDispatch();
  const invoke = useActionInvoker();
  const [sample, setSample] = useState(null);
  const {openModal, dismissModal} = useModal();
  const user = useReduxState((state) => selectUser(state));

  const showSample = () => {
    requestAnswer({
      query: [model_name, '::distinct', attribute.attribute_name]
    })(dispatch).then(({answer}) => setSample(answer));
  };

  useEffect(() => {
    setSample(null);
  }, [attribute]);

  const classes = useStyles();

  const executeAction = useCallback(
    (action) => {
      action
        .then(({models}) => {
          invoke(addTemplatesAndDocuments(models));
        })
        .catch((err) => {
          invoke(showMessages(err));
        });
    },
    [invoke, addTemplatesAndDocuments, showMessages]
  );

  const handleEditAttribute = useCallback(
    (params) => {
      dismissModal();
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

  const isAdminUser = useMemo(() => {
    if (!user || 0 === Object.keys(user).length) return false;

    return isAdmin(user, CONFIG.project_name);
  }, [user, CONFIG.project_name]);

  if (!attribute) return null;

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
          {attribute.attribute_type == 'string' && (
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
