import * as React from 'react'
import Box from '@mui/system/Box'
import IconButton from '@mui/material/IconButton';
import Image from 'next/image';
import arrowUpRightLight from '/public/images/icons/arrow-up-right-light.svg';
import Link from './link';

export default function LinkoutButton({
  link,
  size,
  tooltip
}: {
  link: string;
  size: string;
  tooltip?: string;
}) {
  let iconSize = size == 'small' ? 27 : 39;
  let image = <Image
    style={{ margin: size == 'small' ? '0px' : '4.5px' }}
    width={iconSize}
    height={iconSize}
    src={arrowUpRightLight}
    alt='Arrow pointing up-right'
  />
  return (
    <IconButton sx={{ background: 'black' }}>
      {
        link ?  <Link style={{ display: 'flex' }}
          tooltip={!!tooltip}
          tooltipContent={tooltip}
          target='_blank'
          href={link}>{image}</Link> : image
      }
    </IconButton>
  )
}
