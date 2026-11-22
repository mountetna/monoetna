import React, {useState} from 'react';

import Button from '@material-ui/core/Button';
import MenuItem from '@material-ui/core/MenuItem';
import Tab from '@material-ui/core/Tab';
import Tabs from '@material-ui/core/Tabs';
import TextField from '@material-ui/core/TextField';
import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';

import {json_post} from 'etna-js/utils/fetch';
import {magmaPath} from 'etna-js/api/magma_api';

import {
  TYPES,
  TYPE_HELP,
  countPill,
  emptyModalityForm,
  fieldGroups,
  formBody,
  formError,
  formFromRow,
  modalityFormBody,
  modalityFormFromRow,
  newFieldForm,
  queryText,
  showUnderOptions
} from './field_form';
import {requestError} from './request_error';

const MUTED = 'rgba(0,0,0,0.54)';
const BORDER = '1px solid rgba(0,0,0,0.12)';
const FAINT = '1px solid rgba(0,0,0,0.06)';

const useStyles = makeStyles(() => ({
  root: {flex: 1, minWidth: 0, overflow: 'auto', padding: '8px 16px'},
  tabs: {minHeight: 36, borderBottom: BORDER},
  tab: {minHeight: 36, minWidth: 0, padding: '6px 16px', textTransform: 'none'},
  intro: {margin: '10px 0 4px', maxWidth: 900, fontSize: 13, color: MUTED},
  list: {maxWidth: 900, textAlign: 'left', marginTop: 12},
  groupHeader: {
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'space-between',
    padding: '14px 12px 8px',
    fontSize: 13,
    fontWeight: 600,
    color: 'rgba(0,0,0,0.7)',
    borderBottom: BORDER
  },
  pill: {
    flexShrink: 0,
    padding: '1px 10px',
    borderRadius: 10,
    background: 'rgba(0,0,0,0.06)',
    color: MUTED,
    fontSize: 11,
    fontWeight: 400,
    whiteSpace: 'nowrap'
  },
  row: {
    display: 'flex',
    alignItems: 'center',
    gap: 10,
    padding: '10px 12px',
    fontSize: 13,
    borderBottom: FAINT,
    '&:hover': {background: 'rgba(0,0,0,0.03)'},
    '&:hover $rowActions, &:focus-within $rowActions': {opacity: 1}
  },
  rowText: {flex: 1, minWidth: 0},
  rowLabel: {color: 'rgba(0,0,0,0.87)', overflowWrap: 'anywhere'},
  modalityName: {fontWeight: 600},
  query: {
    marginTop: 1,
    fontFamily: 'monospace',
    fontSize: 11,
    color: MUTED,
    whiteSpace: 'pre-wrap',
    wordBreak: 'break-word'
  },
  rowActions: {opacity: 0, whiteSpace: 'nowrap', transition: 'opacity 0.1s'},
  link: {
    minWidth: 0,
    padding: '0 4px',
    fontSize: 12,
    fontWeight: 400,
    textTransform: 'none',
    color: MUTED,
    '&:hover': {color: 'rgba(0,0,0,0.87)', background: 'transparent', textDecoration: 'underline'}
  },
  none: {padding: '12px', fontSize: 13, fontStyle: 'italic', color: 'rgba(0,0,0,0.38)'},
  queryInput: {fontFamily: 'monospace', fontSize: 12},
  oneLine: {whiteSpace: 'nowrap'},
  form: {display: 'flex', flexDirection: 'column', maxWidth: 640, gap: 12, marginTop: 12},
  error: {marginTop: 8},
  actions: {display: 'flex', gap: 8}
}));

const LinkButton = ({label, onClick}) => {
  const classes = useStyles();
  return <Button size='small' className={classes.link} onClick={onClick}>{label}</Button>;
};

const RowActions = ({canEdit, onEdit, onDelete}) => {
  const classes = useStyles();
  if (!canEdit) return null;

  return (
    <div className={classes.rowActions}>
      <LinkButton label='Edit' onClick={onEdit} />
      <LinkButton label='Delete' onClick={onDelete} />
    </div>
  );
};

const ModalityList = ({modalities, error, canEdit, onAdd, onEdit, onDelete}) => {
  const classes = useStyles();
  const rows = modalities || [];

  return (
    <React.Fragment>
      <Typography className={classes.intro}>A modality is an assay family. A project has it when this model has records.</Typography>
      {canEdit && <Button variant='outlined' color='primary' onClick={onAdd}>+ Add</Button>}
      {error && <Typography color='error' className={classes.error}>{error}</Typography>}
      <div className={classes.list}>
        {rows.length === 0 && <div className={classes.none}>Nothing yet</div>}
        {rows.map((row) => (
          <div className={classes.row} key={row.id}>
            <div className={classes.rowText}>
              <div className={`${classes.rowLabel} ${classes.modalityName}`}>{row.name}</div>
              <div className={classes.query}>{row.model}</div>
            </div>
            <RowActions canEdit={canEdit} onEdit={() => onEdit(row)} onDelete={() => onDelete(row)} />
          </div>
        ))}
      </div>
    </React.Fragment>
  );
};

const FieldList = ({fieldRows, modalities, error, canEdit, onAdd, onEdit, onDelete}) => {
  const classes = useStyles();
  const groups = fieldGroups(fieldRows, modalities);
  const empty = (fieldRows || []).length === 0 && (modalities || []).length === 0;

  return (
    <React.Fragment>
      <Typography className={classes.intro}>A field is a question. Global fields always show. The others show when that modality is selected.</Typography>
      {canEdit && <Button variant='outlined' color='primary' onClick={onAdd}>+ Add</Button>}
      {error && <Typography color='error' className={classes.error}>{error}</Typography>}
      <div className={classes.list}>
        {empty && <div className={classes.none}>Nothing yet</div>}
        {!empty && groups.map((group) => (
          <React.Fragment key={group.key}>
            <div className={classes.groupHeader}>
              <span>{group.title}</span>
              <span className={classes.pill}>{countPill(group.count, 'field')}</span>
            </div>
            {group.rows.map((row) => (
              <div className={classes.row} key={row.id}>
                <div className={classes.rowText}>
                  <div className={classes.rowLabel}>{row.label}</div>
                  <div className={classes.query}>{queryText(row.query)}</div>
                </div>
                <span className={classes.pill}>{row.kind}</span>
                <RowActions canEdit={canEdit} onEdit={() => onEdit(row)} onDelete={() => onDelete(row)} />
              </div>
            ))}
          </React.Fragment>
        ))}
      </div>
    </React.Fragment>
  );
};

const FieldsAdmin = ({fieldRows, modalities, canEdit = false, onChange}) => {
  const classes = useStyles();
  const [tab, setTab] = useState('modalities');
  const [form, setForm] = useState(null);
  const [error, setError] = useState(null);

  const done = () => { setForm(null); setError(null); onChange(); };
  const fail = (failure) => requestError(failure).then(setError);
  const cancel = () => { setForm(null); setError(null); };

  const save = () => {
    if (form.target === 'modality') {
      const body = modalityFormBody(form);
      const request = form.id == null
        ? json_post(magmaPath('cross_project_search/modalities'), body)
        : json_post(magmaPath(`cross_project_search/modalities/${form.id}`), body);
      request.then(done).catch(fail);
      return;
    }
    const syntax = formError(form);
    if (syntax) { setError(syntax); return; }
    const body = formBody(form);
    const request = form.id == null
      ? json_post(magmaPath('cross_project_search/fields'), body)
      : json_post(magmaPath(`cross_project_search/fields/${form.id}`), body);
    request.then(done).catch(fail);
  };

  const removeField = (row) => {
    if (!confirm(`Delete "${row.label}"?`)) return;
    setError(null);
    json_post(magmaPath(`cross_project_search/fields/${row.id}/delete`)).then(onChange).catch(fail);
  };

  const removeModality = (row) => {
    if (!confirm(`Delete "${row.name}"?`)) return;
    setError(null);
    json_post(magmaPath(`cross_project_search/modalities/${row.id}/delete`)).then(onChange).catch(fail);
  };

  if (form && form.target === 'modality') {
    return (
      <div className={classes.root}>
        <Typography variant='h6'>{form.id == null ? 'Add modality' : 'Edit modality'}</Typography>
        <div className={classes.form}>
          <TextField label='Name' value={form.name} onChange={(e) => setForm({...form, name: e.target.value})} />
          <TextField
            label='Model'
            helperText='Template model name.'
            value={form.model}
            onChange={(e) => setForm({...form, model: e.target.value})}
          />
          {error && <Typography color='error' className={classes.error}>{error}</Typography>}
          <div className={classes.actions}>
            <Button variant='contained' color='primary' onClick={save}>Save</Button>
            <Button onClick={cancel}>Cancel</Button>
          </div>
        </div>
      </div>
    );
  }

  if (form) {
    return (
      <div className={classes.root}>
        <Typography variant='h6'>{form.id == null ? 'Add field' : 'Edit field'}</Typography>
        <div className={classes.form}>
          <TextField select label='Type' value={form.kind} onChange={(e) => setForm({...form, kind: e.target.value})}>
            {TYPES.map((type) => <MenuItem key={type.value} value={type.value}>{type.label}</MenuItem>)}
          </TextField>
          <TextField label='Label' value={form.label} onChange={(e) => setForm({...form, label: e.target.value})} />
          <TextField
            label='Query'
            helperText={TYPE_HELP[form.kind]}
            FormHelperTextProps={{className: classes.oneLine}}
            multiline
            minRows={3}
            value={form.queryText}
            onChange={(e) => setForm({...form, queryText: e.target.value})}
            InputProps={{className: classes.queryInput}}
          />
          <TextField select label='Show under' value={form.modalityId} onChange={(e) => setForm({...form, modalityId: e.target.value})}>
            {showUnderOptions(modalities).map((option) => (
              <MenuItem key={option.value || 'global'} value={option.value}>{option.label}</MenuItem>
            ))}
          </TextField>
          <TextField label='Position' type='number' value={form.position} onChange={(e) => setForm({...form, position: e.target.value})} />
          {error && <Typography color='error' className={classes.error}>{error}</Typography>}
          <div className={classes.actions}>
            <Button variant='contained' color='primary' onClick={save}>Save</Button>
            <Button onClick={cancel}>Cancel</Button>
          </div>
        </div>
      </div>
    );
  }

  const showTab = (next) => { setTab(next); setError(null); };

  return (
    <div className={classes.root}>
      <Tabs className={classes.tabs} value={tab} onChange={(e, next) => showTab(next)} indicatorColor='primary' textColor='primary'>
        <Tab className={classes.tab} disableRipple value='modalities' label='Modalities' />
        <Tab className={classes.tab} disableRipple value='fields' label='Fields' />
      </Tabs>
      {tab === 'modalities' ? (
        <ModalityList
          modalities={modalities}
          error={error}
          canEdit={canEdit}
          onAdd={() => { setError(null); setForm({...emptyModalityForm(), target: 'modality'}); }}
          onEdit={(row) => { setError(null); setForm({...modalityFormFromRow(row), target: 'modality'}); }}
          onDelete={removeModality}
        />
      ) : (
        <FieldList
          fieldRows={fieldRows}
          modalities={modalities}
          error={error}
          canEdit={canEdit}
          onAdd={() => { setError(null); setForm({...newFieldForm(null), target: 'field'}); }}
          onEdit={(row) => { setError(null); setForm({...formFromRow(row), target: 'field'}); }}
          onDelete={removeField}
        />
      )}
    </div>
  );
};

export default FieldsAdmin;
