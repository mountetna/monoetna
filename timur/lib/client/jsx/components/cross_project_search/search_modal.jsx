import React, {useEffect, useRef, useState} from 'react';

import Button from '@material-ui/core/Button';
import Checkbox from '@material-ui/core/Checkbox';
import Popover from '@material-ui/core/Popover';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogTitle from '@material-ui/core/DialogTitle';
import IconButton from '@material-ui/core/IconButton';
import LinearProgress from '@material-ui/core/LinearProgress';
import Tab from '@material-ui/core/Tab';
import Tabs from '@material-ui/core/Tabs';
import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';
import ChevronLeftIcon from '@material-ui/icons/ChevronLeft';
import ChevronRightIcon from '@material-ui/icons/ChevronRight';
import ExpandMoreIcon from '@material-ui/icons/ExpandMore';

import {json_get, json_post} from 'etna-js/utils/fetch';
import {downloadBlob} from 'etna-js/utils/blob';
import {magmaPath} from 'etna-js/api/magma_api';

import {
  COUNT,
  PROJECT,
  activeModalities,
  liveModalityOptions,
  modalityGroup,
  MODALITY_ID,
  sidebarFacets,
  buildCsv,
  buildProjects,
  csvRows,
  parseHiddenColumns,
  readHiddenColumns,
  writeHiddenColumns,
  toggleHidden,
  toggleableColumns,
  visibleColumns,
  filterProjects,
  summaryText,
  liveFacetOptions,
  totalsRow,
  tableColumns,
  tableRows,
  toggleValue
} from './facet_filter';
import FieldsAdmin from './fields_admin';
import {normalizeModalities, normalizeFields} from './field_form';
import {requestError} from './request_error';

const useStyles = makeStyles(() => ({
  paper: {fontFamily: 'Open Sans, sans-serif', height: '85vh'},
  title: {display: 'flex', alignItems: 'center', gap: 24, paddingBottom: 0},
  tabs: {minHeight: 36},
  tab: {minHeight: 36, minWidth: 0, padding: '6px 16px', textTransform: 'none'},
  content: {display: 'flex', padding: 0, borderTop: '1px solid rgba(0,0,0,0.12)'},
  sidebar: {
    width: 280,
    flexShrink: 0,
    overflowY: 'auto',
    borderRight: '1px solid rgba(0,0,0,0.12)',
    padding: '8px 16px'
  },
  sidebarCollapsed: {
    width: 40,
    flexShrink: 0,
    borderRight: '1px solid rgba(0,0,0,0.12)',
    padding: '4px 0',
    textAlign: 'center'
  },
  sidebarHeader: {display: 'flex', alignItems: 'center', justifyContent: 'space-between'},
  sidebarTitle: {fontWeight: 600, fontSize: 16},
  clear: {textTransform: 'none', fontSize: 12},
  facet: {borderBottom: '1px solid rgba(0,0,0,0.08)'},
  facetHeader: {
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'space-between',
    width: '100%',
    background: 'none',
    border: 0,
    padding: '10px 0',
    cursor: 'pointer',
    fontFamily: 'inherit',
    fontSize: 14,
    textAlign: 'left'
  },
  facetHeaderActive: {fontWeight: 600},
  options: {paddingBottom: 8, maxHeight: 260, overflowY: 'auto'},
  option: {display: 'flex', alignItems: 'center', fontSize: 14, cursor: 'pointer'},
  checkbox: {padding: 4},
  nested: {paddingLeft: 24},
  nestedLabel: {fontSize: 11, color: 'rgba(0,0,0,0.54)', textTransform: 'uppercase', letterSpacing: 0.4, padding: '4px 0 0 4px'},
  empty: {fontSize: 13, color: 'rgba(0,0,0,0.54)', padding: '0 4px 8px'},
  main: {flex: 1, minWidth: 0, display: 'flex', flexDirection: 'column', padding: '8px 16px'},
  resultsBar: {display: 'flex', alignItems: 'center', justifyContent: 'space-between', marginBottom: 8},
  tableWrap: {overflowX: 'auto', overflowY: 'auto', flex: 1},
  table: {borderCollapse: 'collapse', width: '100%', minWidth: 'max-content', fontSize: 14},
  cell: {
    borderBottom: '1px solid rgba(0,0,0,0.12)',
    padding: '6px 8px',
    textAlign: 'left',
    verticalAlign: 'top',
    wordBreak: 'break-word'
  },
  nowrap: {whiteSpace: 'nowrap', minWidth: 140},
  filterCell: {minWidth: 160, maxWidth: 320},
  countCell: {minWidth: 90, whiteSpace: 'nowrap'},
  head: {position: 'sticky', top: 0, background: '#fff', fontWeight: 600},
  stickyProject: {position: 'sticky', left: 0, background: '#fff', zIndex: 2},
  zero: {color: 'rgba(0,0,0,0.38)'},
  muted: {color: 'rgba(0,0,0,0.38)'},
  totalRow: {fontWeight: 600, background: '#f7f7f7', '& td': {background: '#f7f7f7'}},
  colsPaper: {padding: '8px 12px', minWidth: 220, maxHeight: 400},
  colsHint: {fontSize: 11, color: 'rgba(0,0,0,0.54)', marginTop: 8},
  colsHeader: {display: 'flex', justifyContent: 'space-between', alignItems: 'center'},
  barRight: {display: 'flex', alignItems: 'center', gap: 8},
  error: {margin: 16}
}));

const columnClass = (classes, column) => {
  if (column.kind === PROJECT) return `${classes.nowrap} ${classes.stickyProject}`;
  return column.kind === COUNT ? classes.countCell : classes.filterCell;
};

const Option = ({value, count, checked, onToggle}) => {
  const classes = useStyles();
  return (
    <label className={`${classes.option} ${count === 0 && !checked ? classes.zero : ''}`}>
      <Checkbox
        className={classes.checkbox}
        color='primary'
        size='small'
        checked={checked}
        onChange={onToggle}
      />
      <span>{value} ({count})</span>
    </label>
  );
};

const FacetShell = ({title, tickedCount, open, onToggleOpen, children}) => {
  const classes = useStyles();
  return (
    <div className={classes.facet}>
      <button
        type='button'
        className={`${classes.facetHeader} ${tickedCount ? classes.facetHeaderActive : ''}`}
        aria-expanded={open}
        onClick={onToggleOpen}
      >
        <span>{title}{tickedCount ? ` (${tickedCount})` : ''}</span>
        {open ? <ExpandMoreIcon fontSize='small' /> : <ChevronRightIcon fontSize='small' />}
      </button>
      {open && <div className={classes.options}>{children}</div>}
    </div>
  );
};

const Facet = ({title, options, selected, open, onToggleOpen, onToggleValue}) => {
  const classes = useStyles();
  return (
    <FacetShell title={title} tickedCount={selected.length} open={open} onToggleOpen={onToggleOpen}>
      {options.length === 0 && <div className={classes.empty}>No values</div>}
      {options.map(({value, count}) => (
        <Option
          key={value}
          value={value}
          count={count}
          checked={selected.includes(value)}
          onToggle={() => onToggleValue(value)}
        />
      ))}
    </FacetShell>
  );
};

// Modalities with their nested filter values always listed underneath. A field
// label is shown only when a modality has more than one nested filter.
const ModalityFacet = ({title, options, selected, open, onToggleOpen, onToggleValue}) => {
  const classes = useStyles();
  const ticked = (name) => selected[name] || [];
  const tickedCount = options.reduce(
    (sum, option) => sum + option.groups.reduce((inner, group) => inner + ticked(group.name).length, 0),
    ticked(MODALITY_ID).length
  );
  return (
    <FacetShell title={title} tickedCount={tickedCount} open={open} onToggleOpen={onToggleOpen}>
      {options.length === 0 && <div className={classes.empty}>No values</div>}
      {options.map((option) => (
        <React.Fragment key={option.value}>
          <Option
            value={option.value}
            count={option.count}
            checked={ticked(MODALITY_ID).includes(option.value)}
            onToggle={() => onToggleValue(MODALITY_ID, option.value)}
          />
          {option.groups.filter((group) => group.options.length > 0).map((group) => (
            <div key={group.name} className={classes.nested}>
              {option.groups.length > 1 && <div className={classes.nestedLabel}>{group.label}</div>}
              {group.options.map(({value, count}) => (
                <Option
                  key={value}
                  value={value}
                  count={count}
                  checked={ticked(group.name).includes(value)}
                  onToggle={() => onToggleValue(group.name, value)}
                />
              ))}
            </div>
          ))}
        </React.Fragment>
      ))}
    </FacetShell>
  );
};

const SearchView = ({data, selected, setSelected}) => {
  const classes = useStyles();
  const [expanded, setExpanded] = useState({});
  const [sidebarOpen, setSidebarOpen] = useState(true);
  const [hiddenRaw, setHiddenRaw] = useState(() => readHiddenColumns());
  const [colsAnchor, setColsAnchor] = useState(null);
  const tableWrapRef = useRef(null);

  // Always start with the Project column in view.
  useEffect(() => {
    if (tableWrapRef.current) tableWrapRef.current.scrollLeft = 0;
  }, [data]);

  const fields = (data && data.fields) || [];
  const projects = (data && data.projects) || [];
  const modalities = (data && data.modalities) || [];
  const facets = sidebarFacets(fields, modalities);
  const group = modalityGroup(fields, modalities);

  // A few dozen projects: cheap enough to derive on every render.
  const options = facets.reduce((all, facet) => ({
    ...all,
    [facet.name]: facet.name === MODALITY_ID
      ? liveModalityOptions(modalities, projects, selected, fields)
      : liveFacetOptions(facet.name, projects, selected, group)
  }), {});
  const shown = filterProjects(projects, selected, group);
  const allColumns = tableColumns(fields, modalities, activeModalities(selected, fields, modalities));
  const toggleable = toggleableColumns(allColumns);
  const hidden = parseHiddenColumns(hiddenRaw, fields, modalities);
  const hiddenNow = hidden.filter((id) => allColumns.some((column) => column.id === id));
  const columns = visibleColumns(allColumns, hidden);
  const rows = tableRows(shown, columns);

  const setHidden = (next) => {
    const text = JSON.stringify(next);
    setHiddenRaw(text);
    writeHiddenColumns(text);
  };

  const downloadCsv = () => downloadBlob({
    data: buildCsv(columns, csvRows(shown, columns)),
    filename: 'cross_project_search.csv',
    contentType: 'text/csv;charset=utf-8'
  });

  return (
    <DialogContent className={classes.content}>
      {sidebarOpen ? (
        <div className={classes.sidebar}>
          <div className={classes.sidebarHeader}>
            <Typography className={classes.sidebarTitle}>Filters</Typography>
            <IconButton size='small' aria-label='Collapse filters' onClick={() => setSidebarOpen(false)}>
              <ChevronLeftIcon fontSize='small' />
            </IconButton>
          </div>
          {Object.keys(selected).length > 0 && (
            <Button className={classes.clear} color='primary' size='small' onClick={() => setSelected({})}>
              Clear all filters
            </Button>
          )}
          {facets.map((facet) => {
            const toggleOpen = () => setExpanded((current) => ({...current, [facet.name]: !current[facet.name]}));
            return facet.name === MODALITY_ID ? (
              <ModalityFacet
                key={facet.name}
                title={facet.label}
                options={options[facet.name] || []}
                selected={selected}
                open={Boolean(expanded[facet.name])}
                onToggleOpen={toggleOpen}
                onToggleValue={(name, value) => setSelected((current) => toggleValue(current, name, value))}
              />
            ) : (
              <Facet
                key={facet.name}
                title={facet.label}
                options={options[facet.name] || []}
                selected={selected[facet.name] || []}
                open={Boolean(expanded[facet.name])}
                onToggleOpen={toggleOpen}
                onToggleValue={(value) => setSelected((current) => toggleValue(current, facet.name, value))}
              />
            );
          })}
        </div>
      ) : (
        <div className={classes.sidebarCollapsed}>
          <IconButton size='small' aria-label='Expand filters' onClick={() => setSidebarOpen(true)}>
            <ChevronRightIcon fontSize='small' />
          </IconButton>
        </div>
      )}
      <div className={classes.main}>
        {data && (
          <React.Fragment>
            <div className={classes.resultsBar}>
              <Typography variant='body2'>{summaryText(shown, projects)}</Typography>
              <div className={classes.barRight}>
                <Button variant='outlined' color='primary' onClick={(event) => setColsAnchor(event.currentTarget)}>
                  Columns{hiddenNow.length ? ` (${toggleable.length - hiddenNow.length}/${toggleable.length})` : ''}
                </Button>
                <Button variant='outlined' color='primary' onClick={downloadCsv}>Export CSV</Button>
              </div>
              <Popover
                open={Boolean(colsAnchor)}
                anchorEl={colsAnchor}
                onClose={() => setColsAnchor(null)}
                anchorOrigin={{vertical: 'bottom', horizontal: 'right'}}
                transformOrigin={{vertical: 'top', horizontal: 'right'}}
                classes={{paper: classes.colsPaper}}
              >
                <div className={classes.colsHeader}>
                  <Typography variant='subtitle2'>Columns</Typography>
                  {hiddenNow.length > 0 && <Button size='small' color='primary' className={classes.clear} onClick={() => setHidden(hidden.filter((id) => !hiddenNow.includes(id)))}>Show all</Button>}
                </div>
                {toggleable.map((column) => (
                  <label key={column.id} className={classes.option}>
                    <Checkbox
                      className={classes.checkbox}
                      color='primary'
                      size='small'
                      checked={!hidden.includes(column.id)}
                      onChange={() => setHidden(toggleHidden(hidden, column.id))}
                    />
                    <span>{column.label}</span>
                  </label>
                ))}
                <div className={classes.colsHint}>Order follows Position (edit in Fields view).</div>
              </Popover>
            </div>
            <div className={classes.tableWrap} ref={tableWrapRef}>
              <table className={classes.table}>
                <thead>
                  <tr>
                    {columns.map((column) => (
                      <th key={column.id} className={`${classes.cell} ${classes.head} ${columnClass(classes, column)}`}>{column.label}</th>
                    ))}
                  </tr>
                </thead>
                <tbody>
                  {rows.map((row, index) => (
                    <tr key={shown[index].project_name}>
                      {row.map((text, cellIndex) => (
                        <td
                          key={columns[cellIndex].id}
                          className={`${classes.cell} ${columnClass(classes, columns[cellIndex])} ${text === 'n/a' ? classes.muted : ''}`}
                        >
                          {text}
                        </td>
                      ))}
                    </tr>
                  ))}
                  {shown.length > 0 && (
                    <tr className={classes.totalRow}>
                      {totalsRow(shown, columns).map((text, cellIndex) => (
                        <td key={columns[cellIndex].id} className={`${classes.cell} ${columnClass(classes, columns[cellIndex])}`}>{text}</td>
                      ))}
                    </tr>
                  )}
                </tbody>
              </table>
            </div>
          </React.Fragment>
        )}
      </div>
    </DialogContent>
  );
};

const CrossProjectSearchModal = ({open, canEdit = false, onClose}) => {
  const classes = useStyles();
  const [data, setData] = useState(null);
  const [view, setView] = useState('search');
  const [version, setVersion] = useState(0);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [selected, setSelected] = useState({});

  // The rows are set before the batch runs so the Fields view still works
  // when a query fails.
  useEffect(() => {
    if (!open) return undefined;

    let cancelled = false;
    setLoading(true);
    setError(null);
    Promise.all([
      json_get(magmaPath('cross_project_search/fields')),
      json_get(magmaPath('cross_project_search/modalities'))
    ])
      .then(([{fields: fieldRows}, {modalities: modalityRows}]) => {
        const fields = normalizeFields(fieldRows);
        const modalities = normalizeModalities(modalityRows);
        const base = {fieldRows, modalityRows, fields, modalities, projects: []};
        if (!cancelled) setData(base);
        if (fields.length === 0 && modalities.length === 0) return null;
        const queries = fields.map((field) => field.query).concat(modalities.map((modality) => modality.query));
        return json_post(magmaPath('cross_project_search/query'), {queries}).then(({projects}) => {
          if (!cancelled) setData({...base, projects: buildProjects(fields, modalities, projects)});
        });
      })
      .catch((loadError) => requestError(loadError).then((message) => {
        if (!cancelled) setError(message);
      }))
      .finally(() => {
        if (!cancelled) setLoading(false);
      });

    return () => {
      cancelled = true;
    };
  }, [open, version]);

  return (
    <Dialog classes={{paper: classes.paper}} maxWidth='xl' fullWidth open={open} onClose={onClose}>
      <DialogTitle disableTypography className={classes.title}>
        <Typography variant='h6'>Cross-project search</Typography>
        <Tabs className={classes.tabs} value={view} onChange={(e, next) => setView(next)} indicatorColor='primary' textColor='primary'>
          <Tab className={classes.tab} disableRipple value='search' label='Search' />
          <Tab className={classes.tab} disableRipple value='fields' label='Fields' />
        </Tabs>
      </DialogTitle>
      {loading && <LinearProgress />}
      {error && <Typography color='error' className={classes.error}>{error}</Typography>}
      {view === 'fields' ? (
        <DialogContent className={classes.content}>
          {data && (
            <FieldsAdmin
              fieldRows={data.fieldRows}
              modalities={data.modalityRows}
              canEdit={canEdit}
              onChange={() => setVersion((current) => current + 1)}
            />
          )}
        </DialogContent>
      ) : (
        <SearchView data={data} selected={selected} setSelected={setSelected} />
      )}
      <DialogActions>
        <Button color='secondary' onClick={onClose}>Close</Button>
      </DialogActions>
    </Dialog>
  );
};

export default CrossProjectSearchModal;
