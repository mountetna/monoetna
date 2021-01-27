// Framework libraries.
import * as React from 'react';
import {connect} from 'react-redux';

import Nav from 'etna-js/components/Nav';
import Link from 'etna-js/components/link';
import {selectUser} from 'etna-js/selectors/user-selector';

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
  <div id='vulcan-logo' className='etna-vulcan-logo'>
    <div className='image' />
    <Halo radius={25} />
  </div>
));

const getTabs = () => ({
  workflows: '/workflows',
  help: 'https://mountetna.github.io/vulcan.html'
});

const ModeBar = ({mode, user}) => (
  <div id='nav'>
    {Object.entries(getTabs()).map(([tab_name, route]) => (
      <div
        key={tab_name}
        className={`nav_tab ${mode == tab_name ? 'selected' : ''}`}
      >
        <Link link={route}>{tab_name}</Link>
      </div>
    ))}
  </div>
);

const VulcanNav = ({mode, user}) => (
  <Nav user={user} logo={Logo} app='vulcan'>
    {mode !== 'home' && <ModeBar mode={mode} user={user} />}
  </Nav>
);

export default connect((state) => ({
  user: selectUser(state)
}))(VulcanNav);
