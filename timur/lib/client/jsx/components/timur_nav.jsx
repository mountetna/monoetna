// Framework libraries.
import * as React from 'react';
import {connect} from 'react-redux';

import IconButton from '@material-ui/core/IconButton';
import Tooltip from '@material-ui/core/Tooltip';
import SearchIcon from '@material-ui/icons/Search';
import {makeStyles} from '@material-ui/core/styles';

import Nav from 'etna-js/components/Nav';
import IdentifierSearch from './identifier_search';
import Link from 'etna-js/components/link';
import {selectUser} from 'etna-js/selectors/user-selector';
import {isAdmin, isSuperEditor, isSuperViewer} from 'etna-js/utils/janus';
import {useFeatureFlag} from 'etna-js/hooks/useFeatureFlag';
import CrossProjectSearchModal from './cross_project_search/search_modal';

const Halo = ({radius}) => (
  <div className='halo'>
    <svg>
      <circle r={`${radius * 0.7}px`} cx={`${radius}px`} cy={`${radius}px`} />
      {Array(36)
        .fill()
        .map((_, i) => {
          let deg = i * 10;
          let rad = radius * (i % 2 == 0 ? 1.0 : 0.9);
          let x = (r) => Math.cos((Math.PI * deg) / 180) * r + radius;
          let y = (r) => Math.sin((Math.PI * deg) / 180) * r + radius;
          return (
            <path
              key={i}
              className={i % 2 == 0 ? 'long' : 'short'}
              d={`M ${x(rad)} ${y(rad)} L ${x(radius * 0.7)} ${y(
                radius * 0.7
              )}`}
            />
          );
        })}
    </svg>
  </div>
);

const Logo = connect(({exchanges}) => ({exchanges}))(({exchanges}) => (
  <div
    id='timur-logo'
    className={Object.keys(exchanges).length > 0 ? 'throb' : null}
  >
    <div className='image' />
    <Halo radius={25} />
  </div>
));

const getTabs = ({mode, user, canUseAdvanced}) => {
  let tabs = {
    browse: Routes.browse_path(CONFIG.project_name),
    model: Routes.map_path(CONFIG.project_name),
    search: Routes.search_path(CONFIG.project_name),
    query: Routes.query_path(CONFIG.project_name),
    help: `https://mountetna.github.io/timur.html#${mode}`
  };

  if (isAdmin(user, CONFIG.project_name)) {
    tabs['views'] = Routes.views_path(CONFIG.project_name);
  }

  if (canUseAdvanced) {
    tabs['manifests'] = Routes.manifests_path(CONFIG.project_name);
    tabs['plots'] = Routes.plots_path(CONFIG.project_name);
  }

  return tabs;
};

const ModeBar = ({mode, user}) => {
  const canUseAdvanced = useFeatureFlag('timuradvanced');
  return (
    <div id='nav'>
      {Object.entries(getTabs({mode, user, canUseAdvanced})).map(([tab_name, route]) => (
        <div
          key={tab_name}
          className={`nav_tab ${mode == tab_name ? 'selected' : ''}`}
        >
          <Link link={route}>{tab_name}</Link>
        </div>
      ))}
    </div>
  );
};

const useStyles = makeStyles(() => ({
  search: {
    marginLeft: 15
  },
  spacer: {
    flex: 1
  }
}));

const CrossProjectSearchButton = ({canEdit}) => {
  const classes = useStyles();
  const [open, setOpen] = React.useState(false);

  return (
    <React.Fragment>
      <Tooltip title='Cross-project search'>
        <IconButton className={classes.search} aria-label='Cross-project search' onClick={() => setOpen(true)}>
          <SearchIcon />
        </IconButton>
      </Tooltip>
      <CrossProjectSearchModal open={open} canEdit={canEdit} onClose={() => setOpen(false)} />
    </React.Fragment>
  );
};

const TimurNav = ({mode, user}) => {
  const classes = useStyles();
  const superViewer = Boolean(user.permissions && isSuperViewer(user));
  const superEditor = Boolean(user.permissions && isSuperEditor(user));
  const showModeBar = mode !== 'home';

  return (
    <Nav user={user} logo={Logo} app='timur'>
      {superViewer && mode === 'home' && <CrossProjectSearchButton canEdit={superEditor} />}
      {showModeBar && <IdentifierSearch />}
      {showModeBar && <ModeBar mode={mode} user={user} />}
      {!showModeBar && <div className={classes.spacer} />}
    </Nav>
  );
};

export default connect((state) => ({
  user: selectUser(state)
}))(TimurNav);
