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

import {showMessages} from 'etna-js/actions/message_actions';
import {addTemplatesAndDocuments} from 'etna-js/actions/magma_actions';
import {updateAttribute, removeAttribute} from '../../api/magma_api';
import MapHeading from './map_heading';
import EditAttributeModal from './edit_attribute_modal';

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
  handleRemoveAttribute
}) => {
  const classes = useStyles();
  const {openModal} = useModal();

  const handleConfirmRemove = useCallback(() => {
    if (
      confirm(
        'Removing an attribute is not reversible -- are you sure? If you may want to use the attribute later, you can make it "hidden" so users cannot see it.'
      )
    ) {
      handleRemoveAttribute();
    }
  }, [handleRemoveAttribute]);

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

  const handleEditAttribute = useCallback(
    (params) => {
      dismissModal();
      updateAttribute({
        model_name,
        ...params
      })
        .then(({models}) => {
          invoke(addTemplatesAndDocuments(models));
        })
        .catch((err) => {
          invoke(showMessages(err));
        });
    },
    [invoke, model_name, dismissModal, showMessages, addTemplatesAndDocuments]
  );

  const handleRemoveAttribute = useCallback(() => {
    removeAttribute({
      model_name,
      attribute_name: attribute.name
    })
      .then(({models}) => {
        invoke(addTemplatesAndDocuments(models));
      })
      .catch((err) => {
        invoke(showMessages(err));
      });
  }, [attribute, invoke, addTemplatesAndDocuments, showMessages]);

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
