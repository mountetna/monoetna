import React, {useState, useCallback, useEffect, useMemo} from 'react';

import {selectModels} from 'etna-js/selectors/magma';
import {useReduxState} from 'etna-js/hooks/useReduxState';

import {SNAKE_CASE} from '../../utils/edit_map';
import ModalSelect from './modal_select';
import ModelActionsModal, { ModelModalParams } from './model_actions_modal';
import {ShrinkingLabelTextField} from './shrinking_label_text_field';

import Grid from '@material-ui/core/Grid';
import TextField from '@material-ui/core/TextField';
import SvgIcon from '@material-ui/core/SvgIcon';

const linestyle = {'color': '#000000', 'fill':'#000000', 'strokeMiterlimit': 10};
function SplitIcon(props: any) {
  return (
    <SvgIcon {...props}>
      <path
        style={{...linestyle, 'strokeWidth': 0.995255}}
        d="M 4.1523438 0.75195312 A 0.70006999 0.70006999 0 0 0 3.7128906 0.94726562 L 0.8828125 3.6972656 A 0.70006999 0.70006999 0 0 0 0.8828125 4.7011719 L 3.7128906 7.4511719 A 0.70006999 0.70006999 0 0 0 4.8222656 6.6289062 L 4.0839844 5.1992188 L 7.8417969 5.1992188 L 10.451172 11.003906 L 4.0820312 11.003906 L 4.8222656 9.5703125 A 0.70006999 0.70006999 0 0 0 4.3164062 8.5605469 A 0.70006999 0.70006999 0 0 0 3.7128906 8.7480469 L 0.8828125 11.498047 A 0.70006999 0.70006999 0 0 0 0.8828125 12.501953 L 3.7128906 15.251953 A 0.70006999 0.70006999 0 0 0 4.8222656 14.429688 L 4.0820312 12.996094 L 10.451172 12.996094 L 7.8417969 18.800781 L 4.0839844 18.800781 L 4.8222656 17.371094 A 0.70006999 0.70006999 0 0 0 4.3164062 16.359375 A 0.70006999 0.70006999 0 0 0 3.7128906 16.548828 L 0.8828125 19.298828 A 0.70006999 0.70006999 0 0 0 0.8828125 20.302734 L 3.7128906 23.052734 A 0.70006999 0.70006999 0 0 0 4.8222656 22.228516 L 4.0839844 20.800781 L 9.1171875 20.800781 L 12.625 12.996094 L 19.917969 12.996094 L 19.177734 14.429688 A 0.70006999 0.70006999 0 0 0 20.287109 15.251953 L 23.117188 12.501953 A 0.70006999 0.70006999 0 0 0 23.117188 11.498047 L 20.287109 8.7480469 A 0.70006999 0.70006999 0 0 0 19.382812 8.6875 A 0.70006999 0.70006999 0 0 0 19.177734 9.5703125 L 19.917969 11.003906 L 12.625 11.003906 L 9.1171875 3.1992188 L 4.0839844 3.1992188 L 4.8222656 1.7714844 A 0.70006999 0.70006999 0 0 0 4.3164062 0.75976562 A 0.70006999 0.70006999 0 0 0 4.1523438 0.75195312 z"
      />
    </SvgIcon>
  );
}
function ChildrenIcon(props: any) {
  return (
    <SvgIcon {...props}>
      <path
        style={linestyle}
        d="M 19.748047 8.5527344 A 0.70006999 0.70006999 0 0 0 19.382812 8.6875 A 0.70006999 0.70006999 0 0 0 19.177734 9.5703125 L 19.917969 11.005859 L 4.0820312 11.005859 L 4.8222656 9.5703125 A -0.70006999 0.70006999 0 0 0 4.6171875 8.6875 A -0.70006999 0.70006999 0 0 0 3.7128906 8.7480469 L 0.8828125 11.498047 A -0.70006999 0.70006999 0 0 0 0.8828125 12.501953 L 3.7128906 15.251953 A -0.70006999 0.70006999 0 0 0 4.8222656 14.429688 L 4.0820312 12.994141 L 19.917969 12.994141 L 19.177734 14.429688 A 0.70006999 0.70006999 0 0 0 20.287109 15.251953 L 23.117188 12.501953 A 0.70006999 0.70006999 0 0 0 23.117188 11.498047 L 20.287109 8.7480469 A 0.70006999 0.70006999 0 0 0 19.748047 8.5527344 z"
      />
      <path
        style={linestyle}
        d="M 4.1523438 0.75195312 A 0.70006999 0.70006999 0 0 0 3.7128906 0.94726562 L 0.8828125 3.6972656 A 0.70006999 0.70006999 0 0 0 0.8828125 4.7011719 L 3.7128906 7.4511719 A 0.70006999 0.70006999 0 0 0 4.8222656 6.6289062 L 4.078125 5.1875 L 19.921875 5.1875 L 19.177734 6.6289062 A -0.70006999 0.70006999 0 0 0 20.287109 7.4511719 L 23.117188 4.7011719 A -0.70006999 0.70006999 0 0 0 23.117188 3.6972656 L 20.287109 0.94726562 A -0.70006999 0.70006999 0 0 0 19.683594 0.75976562 A -0.70006999 0.70006999 0 0 0 19.177734 1.7714844 L 19.916016 3.1992188 L 4.0839844 3.1992188 L 4.8222656 1.7714844 A 0.70006999 0.70006999 0 0 0 4.3164062 0.75976562 A 0.70006999 0.70006999 0 0 0 4.1523438 0.75195312 z"
      />
      <path
        style={linestyle}
        d="M 4.1523438 16.351562 A 0.70006999 0.70006999 0 0 0 3.7128906 16.548828 L 0.8828125 19.298828 A 0.70006999 0.70006999 0 0 0 0.8828125 20.302734 L 3.7128906 23.052734 A 0.70006999 0.70006999 0 0 0 4.8222656 22.228516 L 4.0839844 20.800781 L 19.916016 20.800781 L 19.177734 22.228516 A -0.70006999 0.70006999 0 0 0 20.287109 23.052734 L 23.117188 20.302734 A -0.70006999 0.70006999 0 0 0 23.117188 19.298828 L 20.287109 16.548828 A -0.70006999 0.70006999 0 0 0 19.683594 16.359375 A -0.70006999 0.70006999 0 0 0 19.177734 17.371094 L 19.921875 18.8125 L 4.078125 18.8125 L 4.8222656 17.371094 A 0.70006999 0.70006999 0 0 0 4.3164062 16.359375 A 0.70006999 0.70006999 0 0 0 4.1523438 16.351562 z"
      />
    </SvgIcon>
  );
}

export default function AddLinkModal({modelName,onSave,open,onClose}: ModelModalParams & {modelName: string}) {
  const [linkAttributeName, setLinkAttributeName] = useState('');
  const [reciprocalModelName, setReciprocalModelName] = useState('');
  const [reciprocalAttributeName, setReciprocalAttributeName] = useState('');
  const [reciprocalLinkType, setReciprocalLinkType] = useState('collection');

  const models = useReduxState((state: any) => selectModels(state));

  const handleOnSave = useCallback(() => {
    onSave({
      linkAttributeName,
      reciprocalModelName,
      reciprocalAttributeName,
      reciprocalLinkType
    });
  }, [
    linkAttributeName,
    reciprocalModelName,
    reciprocalAttributeName,
    reciprocalLinkType
  ]);

  const disabled = !(linkAttributeName &&
      reciprocalModelName &&
      reciprocalAttributeName &&
      reciprocalLinkType) || (
        linkAttributeName==reciprocalAttributeName &&
        modelName==reciprocalModelName
      );

  const reset = useCallback(() => {
    setLinkAttributeName('');
    setReciprocalModelName('');
    setReciprocalAttributeName('');
    setReciprocalLinkType('collection');
  }, []);

  const handleOnCancel = useCallback(() => {
    onClose();
    reset();
  }, []);

  const handleSetReciprocalModelName = useCallback((name) => {
    if (modelName!=name &&
      ((!linkAttributeName && !reciprocalAttributeName) ||
      (linkAttributeName==reciprocalModelName && reciprocalAttributeName==modelName))
    ) {
      setReciprocalAttributeName(modelName)
      setLinkAttributeName(name)
    }
    setReciprocalModelName(name);
  }, [])

  const reciprocalModelNameOptions = useMemo(() => {
    return Object.keys(models);
  }, [models]);

  const reciprocalLinkTypeOptions = ['child', 'collection'];

  const semanticDescription = 
    <div style={{display: 'flex', justifyContent: 'center', alignItems: 'center', gap: '8px', paddingTop: '40px'}}>
      { modelName!=reciprocalModelName ?
        (reciprocalLinkType=='collection' ?
          `Many '${modelName}' records linkable to same '${reciprocalModelName || 'end_model'}' records` :
          `Single '${modelName}' record linkable to each '${reciprocalModelName || 'end_model'}' record`) :
        reciprocalLinkType=='collection' ?
          `Many '${modelName}.${linkAttributeName || 'start_attribute'}' linkable to same '${reciprocalModelName || 'end_model'}.${reciprocalAttributeName || 'end_attribute'}'` :
          `Single '${modelName}.${linkAttributeName || 'start_attribute'}' linkable to each '${reciprocalModelName || 'end_model'}.${reciprocalAttributeName || 'end_attribute'}'`
      }
    </div>

  return (
    <ModelActionsModal onClose={handleOnCancel} open={open} onSave={handleOnSave} title='Add Link' saveDisabled={disabled}>
      <Grid container direction="row" style={{flexWrap: 'nowrap', justifyContent: 'space-between'}}>
        <Grid style={{flex: '1 1 auto', maxWidth: 380, display: 'flex', flexDirection: 'column', gap: '10px'}}>
          <TextField
            disabled
            fullWidth
            id='start-model-name'
            value={modelName}
            label='Start Model'
            InputLabelProps={{shrink: true}}
          />
          <TextField
            disabled
            fullWidth
            id='start-model-type'
            value='link'
            label='Attribute Type'
            InputLabelProps={{shrink: true}}
          />
          <ShrinkingLabelTextField
            id='start-attribute-name'
            label='Attribute Name (snake_case)'
            value={linkAttributeName}
            onChange={(e: React.ChangeEvent<any>) =>
              setLinkAttributeName(e.target.value)
            }
            pattern={SNAKE_CASE}
          />
        </Grid>
        <Grid style={{display: 'flex', width: '40px', justifyContent: 'center', alignItems: 'center'}}>
          {reciprocalLinkType=='collection' ? <SplitIcon fontSize='large'/>: <ChildrenIcon fontSize='large'/>}
        </Grid>
        <Grid style={{flex: '1 1 auto', maxWidth: 380, display: 'flex', flexDirection: 'column', gap: '10px'}}>
          <ModalSelect
            id='reciprocal-model-name'
            value={reciprocalModelName}
            label='End Model'
            onChange={handleSetReciprocalModelName}
            options={reciprocalModelNameOptions}
          />
          <ModalSelect
            id='reciprocal-link-type'
            value={reciprocalLinkType}
            label='Attribute Type'
            onChange={setReciprocalLinkType}
            options={reciprocalLinkTypeOptions}
          />
          <ShrinkingLabelTextField
            id='reciprocal-attribute-name'
            label='Attribute Name (snake_case)'
            value={reciprocalAttributeName}
            onChange={(e: React.ChangeEvent<any>) =>
              setReciprocalAttributeName(e.target.value)
            }
            pattern={SNAKE_CASE}
          />
        </Grid>
      </Grid>
      {semanticDescription}
    </ModelActionsModal>
  );
}
