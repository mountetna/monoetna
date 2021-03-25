import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';

import Icon from 'etna-js/components/icon';
import {STATUS} from '../../../models/steps';
import AnimatedClock from './animated_clock';

export default function StepName({step, status}) {
  let {calculating} = useContext(VulcanContext);

  const icons = {};
  icons[STATUS.COMPLETE] = {
    icon: 'check',
    className: 'light green'
  };
  icons[STATUS.PENDING] = {icon: 'clock', className: 'light'};
  icons[STATUS.ERROR] = {icon: 'times-circle', className: 'light red'};

  let stepStatus = status || STATUS.PENDING;
  let icon = icons[stepStatus];

  let className = `step-status-icon ${icon.className}`;
  let IconComponent = <Icon title={ step.label || step.name } className={className} icon={icon.icon}></Icon>;

  // If the icon is PENDING and also the app state `calculating` == true,
  //   we'll replace the icon with an animated one!
  if (STATUS.PENDING === stepStatus && calculating) {
    IconComponent = <AnimatedClock />;
  }

  return (
    <div className='step-name'>
      <div className='step-status-icon-wrapper'>{IconComponent}</div>
      <div className='step-button'>{step.label || step.name}</div>
    </div>
  );
}
